import pandas as pd
import os
import logging
import argparse
import subprocess
import re
import colorsys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
MAX_LINE_WIDTH = 80
QUALIFIER_START = 21
FEATURE_WIDTH = 15

def convert_dna_to_genbank(dna_file_path, temp_file_prefix="temp_genbank", timeout=300):
    """SnapGene .dna 파일을 .gbk로 변환"""
    logging.info(f"Converting {dna_file_path} to GenBank format")
    snapgene_path = "/Applications/SnapGene.app/Contents/MacOS/SnapGene"
    temp_gb_file = os.path.abspath(f"{temp_file_prefix}.gbk")
    
    if os.path.exists(temp_gb_file):
        logging.info(f"{temp_gb_file} already exists, using existing file")
        return temp_gb_file
    
    cmd = [snapgene_path, "--convert", "GenBank - SnapGene", "--input", dna_file_path, "--output", temp_gb_file]
    logging.debug(f"Executing command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, timeout=timeout, capture_output=True, text=True)
        logging.debug(f"SnapGene output: {result.stdout}")
        if not os.path.exists(temp_gb_file):
            raise FileNotFoundError(f"Output file {temp_gb_file} not generated")
        logging.info(f"Converted {dna_file_path} to {temp_gb_file}")
    except subprocess.TimeoutExpired as e:
        logging.error(f"Conversion timed out after {timeout} seconds: {e.stdout}\n{e.stderr}")
        if os.path.exists(temp_gb_file):
            os.remove(temp_gb_file)
        raise
    except subprocess.CalledProcessError as e:
        logging.error(f"SnapGene CLI failed with exit code {e.returncode}: {e.stderr}")
        if os.path.exists(temp_gb_file):
            os.remove(temp_gb_file)
        raise
    except Exception as e:
        logging.error(f"Error during conversion: {str(e)}")
        if os.path.exists(temp_gb_file):
            os.remove(temp_gb_file)
        raise
    return temp_gb_file

def format_long_qualifier(value, first_line_prefix, subsequent_prefix):
    """긴 qualifier 값을 GenBank 형식에 맞게 포맷팅"""
    lines = []
    remaining_value = value.strip()
    first_line_done = False
    
    while remaining_value:
        if not first_line_done:
            prefix = first_line_prefix
            max_width = MAX_LINE_WIDTH - len(prefix)
            split_value = remaining_value[:max_width]
            lines.append(f"{prefix}\"{split_value}\"")
            remaining_value = remaining_value[max_width:]
            first_line_done = True
        else:
            prefix = subsequent_prefix
            max_width = MAX_LINE_WIDTH - len(prefix)
            split_value = remaining_value[:max_width]
            lines.append(f"{prefix}{split_value}\"")
            remaining_value = remaining_value[max_width:]
        if not remaining_value:
            break
    
    return lines

def hsv_to_hex(h, s, v):
    """HSV 색상을 HEX 코드로 변환"""
    rgb = colorsys.hsv_to_rgb(h / 360.0, s / 100.0, v / 100.0)
    r, g, b = [int(x * 255) for x in rgb]
    return f"#{r:02x}{g:02x}{b:02x}"

def get_color_for_group(group, stage_names):
    """그룹에 따라 채도가 점차 옅어지는 색상 반환"""
    base_hue = 300
    base_saturation = 100
    base_value = 100
    
    if group == "Unknown":
        return "#FF00FF"
    try:
        group_index = stage_names.index(group)
        saturation = max(20, 100 - group_index * 20)
    except ValueError:
        saturation = 100
    
    return hsv_to_hex(base_hue, saturation, base_value)

def parse_mutation_range(position, mutation=None, annotation=None, sheet_name="Mutations", end_position=None):
    """Mutation 열 또는 start/end에서 범위 파싱"""
    pos_str = str(position).replace(',', '').replace('–', '-').rstrip('-')
    ranges = []
    
    if sheet_name == "Mutations":
        start = int(pos_str)
        end = start
        if mutation and pd.notna(mutation) and mutation != "Unknown":
            if 'Δ' in mutation:
                delta = int(mutation.split('Δ')[1].split()[0].replace(',', ''))
                end = start + delta - 1
            elif '→' in mutation and '(' in mutation:
                old, new = mutation.split('→')
                old_count = int(re.search(r'\d+', old).group()) if re.search(r'\d+', old) else 1
                new_count = int(re.search(r'\d+', new).group()) if re.search(r'\d+', new) else 1
                if new_count > old_count:
                    end = start + (new_count - old_count) - 1
                elif new_count < old_count:
                    end = start + (old_count - new_count) - 1
            elif '+' in mutation or '-' in mutation:
                change = mutation.strip('+-')
                end = start + len(change) - 1 if '+' in mutation else start
        start, end = min(start, end), max(start, end)
        ranges.append((start, end))
    
    elif sheet_name == "Missing":
        # start 파싱
        start_parts = pos_str.split('-')
        if len(start_parts) == 2:
            start1, start2 = map(int, start_parts)
        else:
            start1 = start2 = int(pos_str)
        
        # end 파싱
        if end_position and pd.notna(end_position):
            end_str = str(end_position).replace(',', '').replace('–', '-').rstrip('-')
            end_parts = end_str.split('-')
            if len(end_parts) == 2:  # WWW-ZZZ 형식
                end1, end2 = map(int, end_parts)
            else:  # WWW 형식 (단일 값)
                end1 = end2 = int(end_str)
        else:
            end1 = end2 = start2
            end_parts = [str(end1)]
        
        # start와 end가 모두 범위일 경우
        if len(start_parts) == 2 and len(end_parts) == 2:
            ranges.append((start1, end1))  # 첫 번째 시작점과 첫 번째 끝점
            ranges.append((start2, end2))  # 두 번째 시작점과 두 번째 끝점
        # start가 범위이고 end가 단일 값일 경우
        elif len(start_parts) == 2 and len(end_parts) == 1:
            ranges.append((start1, end1))
            ranges.append((start2, end1))
        # start가 단일 값이고 end가 범위일 경우
        elif len(start_parts) == 1 and len(end_parts) == 2:
            ranges.append((start1, end1))
            ranges.append((start1, end2))
        # start와 end가 모두 단일 값일 경우
        else:
            ranges.append((start1, end1))
        
        # 디버깅 로그 추가
        logging.debug(f"Parsed ranges for {sheet_name}: start={pos_str}, end={end_position}, ranges={ranges}")
    
    return ranges

def create_annotated_genbank(dna_file, xlsx_file, output_genbank_file, sheets=["Mutations", "Missing"]):
    """DNA 파일의 기존 GenBank 정보를 유지하며 mutation 정보를 추가"""
    logging.info(f"Processing DNA file: {dna_file}")
    if dna_file.endswith('.dna'):
        temp_gbk = convert_dna_to_genbank(dna_file)
    else:
        temp_gbk = dna_file
    
    with open(temp_gbk, 'r') as f:
        genbank_lines = f.readlines()
    
    feature_start = next(i for i, line in enumerate(genbank_lines) if line.startswith("FEATURES"))
    origin_start = next(i for i, line in enumerate(genbank_lines) if line.startswith("ORIGIN"))
    header_lines = genbank_lines[:feature_start + 1]
    feature_lines = genbank_lines[feature_start + 1:origin_start]
    origin_lines = genbank_lines[origin_start:]
    
    new_feature_lines = []
    for sheet_name in sheets:
        logging.info(f"Reading data from {xlsx_file}, sheet: {sheet_name}")
        df = pd.read_excel(xlsx_file, sheet_name=sheet_name)
        logging.info(f"Loaded {sheet_name} with {len(df)} rows")
        logging.debug(f"Column names in {sheet_name}: {list(df.columns)}")
        
        stage_names = sorted(set(col.split('.')[1].split('/')[0] for col in df.columns if '.' in col and col != 'position' and not col.startswith('Unnamed')))
        if not stage_names:
            logging.warning(f"No stage names detected in {sheet_name}. Using column names as stages.")
            stage_names = [col for col in df.columns if col != 'position' and not col.startswith('Unnamed')]
        logging.info(f"Detected stages in {sheet_name}: {stage_names}")
        
        def find_first_detection(row):
            for stage in stage_names:
                mutation_col = f"mutation.{stage}/output/index.html" if sheet_name == "Mutations" else None
                end_col = f"end.{stage}/output/index.html" if sheet_name == "Missing" else None
                if mutation_col and mutation_col in df.columns and pd.notna(row[mutation_col]):
                    return stage, row[mutation_col]
                if end_col and end_col in df.columns and pd.notna(row[end_col]):
                    return stage, "Missing Coverage"
                for col in df.columns:
                    if stage in col and pd.notna(row[col]):
                        return stage, row.get(mutation_col, "Missing Coverage") if mutation_col else "Missing Coverage"
            return "Unknown", "Unknown"
        
        df[['first_detection', 'first_mutation']] = df.apply(lambda row: pd.Series(find_first_detection(row)), axis=1)
        
        position_col = 'position' if 'position' in df.columns else 'start' if 'start' in df.columns else None
        if not position_col:
            logging.error(f"No 'position' or 'start' column found in {sheet_name}. Skipping.")
            continue
        
        for _, row in df.dropna(subset=[position_col]).iterrows():
            pos_str = str(row[position_col])
            mutation = row['first_mutation']
            stage = row['first_detection']
            annotation = str(row.get('annotation', 'Unknown')) if 'annotation' in df.columns else 'Unknown'
            
            # 동적으로 end 열 찾기
            end_col = f"end.{stage}/output/index.html" if sheet_name == "Missing" else 'end'
            end_position = row.get(end_col, None) if end_col in df.columns else None
            
            if pd.isna(mutation) or mutation.lower() == 'nan':
                mutation = "Missing Coverage"
            if pd.isna(annotation) or annotation.lower() == 'nan':
                annotation = "Unknown"
            
            try:
                ranges = parse_mutation_range(pos_str, mutation, annotation, sheet_name, end_position)
            except (ValueError, IndexError) as e:
                logging.warning(f"Skipping invalid position/mutation in {sheet_name}: {pos_str}/{mutation} - {str(e)}")
                continue
            
            label_prefix = "[Mutation]" if sheet_name == "Mutations" else "[Missing]"
            color = get_color_for_group(stage, stage_names)
            
            for i, (start, end) in enumerate(ranges):
                if sheet_name == "Missing":
                    delta = end - start + 1  # 구간 크기
                    label = f"{label_prefix} Δ{delta} {annotation}".strip()
                    logging.debug(f"Range {i+1} in {sheet_name}: {start}-{end}, delta={delta}, label={label}")
                else:
                    label = f"{label_prefix} {mutation} {annotation}".strip()
                    logging.debug(f"Range {i+1} in {sheet_name}: {start}-{end}, label={label}")
                
                feature_line = f"     misc_feature    {start}..{end}\n"
                new_feature_lines.append(feature_line)
                
                qualifiers = [
                    ("/label", label),
                    ("/note", f"First detected in {stage}"),
                    ("/note", f"color: {color}")
                ]
                for key, value in qualifiers:
                    prefix = f"{' ' * QUALIFIER_START}{key}="
                    subsequent_prefix = " " * QUALIFIER_START
                    formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix)
                    new_feature_lines.extend([line + "\n" for line in formatted_lines])
    
    all_feature_lines = feature_lines + new_feature_lines
    genbank_lines = header_lines + all_feature_lines + origin_lines
    
    with open(output_genbank_file, 'w') as f:
        f.writelines(genbank_lines)
    logging.info(f"Annotated GenBank file saved to: {output_genbank_file}")
    
    if dna_file.endswith('.dna'):
        os.remove(temp_gbk)

def main():
    parser = argparse.ArgumentParser(description="Annotate mutations in GenBank format from DNA and Excel data.")
    parser.add_argument("--dna_file", required=True, help="Path to the DNA file (.dna or .gbk)")
    parser.add_argument("--mutation_xlsx", default="merged_mutant_results.xlsx", help="Path to the Excel file with mutation data")
    parser.add_argument("--output_dir", default="deletion_results", help="Directory to save output files")
    
    args = parser.parse_args()
    
    args.dna_file = os.path.abspath(args.dna_file)
    args.mutation_xlsx = os.path.abspath(args.mutation_xlsx)
    args.output_dir = os.path.abspath(args.output_dir)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    mutation_output_file = os.path.join(args.output_dir, "annotated_dna_with_mutations.gbk")
    create_annotated_genbank(args.dna_file, args.mutation_xlsx, mutation_output_file)

if __name__ == "__main__":
    main()