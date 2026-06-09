import gc
import ctypes
import time
import sys
import platform

def nuclear_ram_flush():
    """Force system to clear old data traces in RAM."""
    gc.collect()
    if platform.system() == 'Linux':
        try:
            # Force Linux to return RAM directly from C-heap to OS
            ctypes.CDLL('libc.so.6').malloc_trim(0)
        except Exception:
            pass
    time.sleep(0.1) # Brief pause for OS to update Page Table

class DualLogger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')
    def write(self, message):
        try:
            self.terminal.write(message)
        except UnicodeEncodeError:
            self.terminal.write(message.encode(self.terminal.encoding, errors='replace').decode(self.terminal.encoding))
        self.log.write(message)
        self.log.flush()
    def flush(self):
        self.terminal.flush()
        self.log.flush()

IUPAC_MAP = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

def generate_iupac_consensus(sequences):
    """Generate degenerate IUPAC consensus sequence from a list of aligned gene segments."""
    if not sequences:
        return ""
    length = len(sequences[0])
    consensus = []
    for i in range(length):
        bases = set()
        for seq in sequences:
            if i < len(seq):
                b = seq[i].upper()
                if b in ['A', 'C', 'G', 'T']:
                    bases.add(b)
        if not bases:
            consensus.append('N')
        else:
            consensus.append(IUPAC_MAP.get(frozenset(bases), 'N'))
    return "".join(consensus)

def generate_batch_analytical_summary(path_master_stats, current_params, language="vi") -> str:
    """
    Automatically analyze primer batches, calculate thermodynamics (Tm, GC, Tm Delta)
    and package into a highly concise analytical summary for the AI Expert.
    """
    import os
    import pandas as pd
    from Bio.SeqUtils import MeltingTemp as mt
    
    if not os.path.exists(path_master_stats):
        return "Không tìm thấy dữ liệu ứng viên để phân tích." if language == "vi" else "No candidate data found for analysis."
        
    try:
        df = pd.read_csv(path_master_stats)
        if df.empty:
            return "Danh sách ứng viên trống." if language == "vi" else "Candidate list is empty."
            
        def calc_tm(seq):
            if pd.isna(seq) or not seq: return 0.0
            # Remove degenerate IUPAC characters if present
            clean_seq = "".join([c for c in str(seq).upper() if c in ['A', 'C', 'G', 'T']])
            if len(clean_seq) < 8: return 0.0
            try:
                # IDT OligoAnalyzer / ThermoFisher-compatible conditions:
                # Na=50mM, Mg=3.0mM, dNTPs=0.8mM, [Oligo]=250nM
                # SantaLucia 1998 (DNA_NN3)
                return round(mt.Tm_NN(
                    clean_seq,
                    nn_table=mt.DNA_NN3,
                    Na=50,
                    Mg=3.0,
                    dNTPs=0.8,
                    dnac1=250,
                    dnac2=0
                ), 1)
            except Exception: return 0.0
            
        def calc_gc(seq):
            if pd.isna(seq) or not seq: return 0.0
            seq_str = str(seq).upper()
            g_c = seq_str.count('G') + seq_str.count('C')
            return round((g_c / len(seq_str)) * 100, 1) if len(seq_str) > 0 else 0.0
            
        def calc_end3_gc(seq):
            """Check if 3' end (last 5 bases) has GC clamp (1-2 G/C is ideal)."""
            if pd.isna(seq) or not seq: return 0
            end = str(seq).upper()[-5:]
            return end.count('G') + end.count('C')
            
        df['Fwd_Tm'] = df['Forward Primer'].apply(calc_tm)
        df['Rev_Tm'] = df['Reverse Primer'].apply(calc_tm)
        df['Tm_Delta'] = (df['Fwd_Tm'] - df['Rev_Tm']).abs()
        df['Fwd_GC'] = df['Forward Primer'].apply(calc_gc)
        df['Rev_GC'] = df['Reverse Primer'].apply(calc_gc)
        df['Fwd_3end_GC'] = df['Forward Primer'].apply(calc_end3_gc)
        df['Rev_3end_GC'] = df['Reverse Primer'].apply(calc_end3_gc)
        
        def get_qc_flag(row):
            flags = []
            # Tm range check (IDT/ThermoFisher-aligned: 55-68°C is optimal for qPCR)
            if row['Fwd_Tm'] < 55.0 or row['Fwd_Tm'] > 68.0: flags.append("[QC_WARN_FWD_TM_OOR]")
            if row['Rev_Tm'] < 55.0 or row['Rev_Tm'] > 68.0: flags.append("[QC_WARN_REV_TM_OOR]")
            if row['Tm_Delta'] > 3.0: flags.append("[QC_WARN_HIGH_DELTA]")
            if row['Fwd_GC'] < 40 or row['Fwd_GC'] > 60: flags.append("[QC_WARN_FWD_GC]")
            if row['Rev_GC'] < 40 or row['Rev_GC'] > 60: flags.append("[QC_WARN_REV_GC]")
            # 3' end GC clamp check (1-2 G/C in last 5 bases)
            fwd_3gc = row.get('Fwd_3end_GC', 0)
            rev_3gc = row.get('Rev_3end_GC', 0)
            if fwd_3gc == 0: flags.append("[QC_WARN_FWD_NO_3GC]")
            if fwd_3gc >= 4: flags.append("[QC_WARN_FWD_3GCC_RICH]")
            if rev_3gc == 0: flags.append("[QC_WARN_REV_NO_3GC]")
            if rev_3gc >= 4: flags.append("[QC_WARN_REV_3GCC_RICH]")
            if not flags: flags.append("[QC_PASS]")
            return " ".join(flags)
            
        df['QC_Flags'] = df.apply(get_qc_flag, axis=1)
        # Save dataframe back with QC flags for downstream traceability
        df.to_csv(path_master_stats, index=False)
        
        min_sens = current_params.get('min_sensitivity', 90.0)
        qualified = df[df['Sensitivity_Percent'] >= min_sens]
        
        if language == "en":
            summary = f"=== BATCH ANALYTICAL PRE-COMPUTATION REPORT ===\n"
            summary += f"Total candidates in batch: {len(df)}\n"
            summary += f"Candidates meeting min sensitivity ({min_sens}%): {len(qualified)}\n\n"
            summary += "--- TOP 15 CANDIDATES RANKED AND ANALYZED (Ranked by Sensitivity & Tm Delta) ---\n"
        else:
            summary = f"=== BÁO CÁO PHÂN TÍCH LÔ MỒI TỰ ĐỘNG (BATCH ANALYTICAL PRE-COMPUTATION REPORT) ===\n"
            summary += f"Tổ̉ng số ứng viên trong lô: {len(df)}\n"
            summary += f"Số ứng viên đạt ngưỡng độ nhạy tối thiểu ({min_sens}%): {len(qualified)}\n\n"
            summary += "--- DANH SÁCH TOP 15 ỨNG VIÊN ĐƯỢC XẾP HẠNG VÀ PHÂN TÍCH (Ranked by Sensitivity & Tm Delta) ---\n"
        
        # Sort candidates prioritizing best sensitivity, then smallest Tm Delta
        df_sorted = df.sort_values(by=['Sensitivity_Percent', 'Tm_Delta'], ascending=[False, True])
        top_candidates = df_sorted.head(15)  # Get top 15 best candidates for cross comparison
        
        for idx, (_, row) in enumerate(top_candidates.iterrows()):
            if language == "en":
                summary += f"Candidate #{idx+1} [Set_ID: {row['Primer Pair']}]:\n"
                summary += f"   - Sensitivity: {row['Sensitivity_Percent']}% | Max Copy: {row['Max_Copy_Number']}\n"
                summary += f"   - Forward Primer: {row['Forward Primer']} (Tm: {row['Fwd_Tm']}°C, GC: {row['Fwd_GC']}%)\n"
                summary += f"   - Reverse Primer: {row['Reverse Primer']} (Tm: {row['Rev_Tm']}°C, GC: {row['Rev_GC']}%)\n"
                summary += f"   - Tm Difference (Tm Delta): {row['Tm_Delta']}°C\n"
            else:
                summary += f"Ứng viên #{idx+1} [Set_ID: {row['Primer Pair']}]:\n"
                summary += f"   - Độ nhạy (Sensitivity): {row['Sensitivity_Percent']}% | Bản sao tối đa (Max Copy): {row['Max_Copy_Number']}\n"
                summary += f"   - Forward Primer: {row['Forward Primer']} (Tm: {row['Fwd_Tm']}°C, GC: {row['Fwd_GC']}%)\n"
                summary += f"   - Reverse Primer: {row['Reverse Primer']} (Tm: {row['Rev_Tm']}°C, GC: {row['Rev_GC']}%)\n"
                summary += f"   - Chênh lệch Tm (Tm Delta): {row['Tm_Delta']}°C\n"
            summary += f"   - QC Flags: {row['QC_Flags']}\n"
            
            # Scan for structural risks
            risks = []
            if row['Tm_Delta'] > 3.0:
                risks.append("Chênh lệch Tm lớn (>3°C)" if language == "vi" else "Large Tm difference (>3°C)")
            if row['Fwd_GC'] < 40 or row['Fwd_GC'] > 60:
                risks.append("GC mồi F ngoài khoảng tối ưu 40-60%" if language == "vi" else "Fwd primer GC out of optimal 40-60% range")
            if row['Rev_GC'] < 40 or row['Rev_GC'] > 60:
                risks.append("GC mồi R ngoài khoảng tối ưu 40-60%" if language == "vi" else "Rev primer GC out of optimal 40-60% range")
                
            if risks:
                if language == "en":
                    summary += f"   - Physical warning: ⚠️ {', '.join(risks)}\n"
                else:
                    summary += f"   - Cảnh báo vật lý: ⚠️ {', '.join(risks)}\n"
            else:
                if language == "en":
                    summary += f"   - Physical warning: ✅ No significant physical risks detected (High stability)\n"
                else:
                    summary += f"   - Cảnh báo vật lý: ✅ Không phát hiện rủi ro vật lý lớn (Độ ổn định cao)\n"
            summary += "\n"
            
        return summary
    except Exception as e:
        return f"Lỗi trong quá trình phân tích số liệu: {e}" if language == "vi" else f"Error in data analysis: {e}"
