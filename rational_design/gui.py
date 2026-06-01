import streamlit as st
import os
import json
import subprocess
import sys
import platform

# --- 1. TRANSLATION SYSTEM (i18n) ---
T = {
    "vi": {
        "title": "🤖 V-Extreme Chuyên gia AI (AI Expert)",
        "subtitle": "Hệ thống thiết kế mồi tối ưu tích hợp AI",
        "welcome": (
            "👋 **Chào bạn! Tôi là V-Extreme Chuyên gia AI (AI Expert).**\n\n"
            "Tôi ở đây để giúp bạn thiết kế mồi PCR một cách dễ dàng nhất. Bạn có thể yêu cầu tôi:\n\n"
            "1-Thiết kế primer mới với dữ liệu từ NCBI\n\n"
            "2-Thiết kế primers mới với dữ liệu local từ thư mục của bạn.\n\n"
            "3-Kiểm tra độ nhạy, độ đặc hiệu của một bộ mồi đã biết\n\n"
            "4-Phân tích kết quả insilico, cho các kết quả đã thiết kế trước đó\n\n"
            "5-Tư vấn tối ưu hoá quy trình wetlab dựa trên các cặp primers/probes đã thiết kế."
        ),
        "lang_selector": "🌐 Chọn ngôn ngữ / Language",
        "reset_search": "🗑️ Reset các ô tìm kiếm",
        "essential_settings": "🔬 Cài Đặt Cơ Bản",
        "ncbi_email": "NCBI Email (Bắt buộc để tải Gen)",
        "save_email": "💾 Lưu email cho lần sau",
        "min_sensitivity": "Độ nhạy tối thiểu (%)",
        "min_sensitivity_help": "Chỉ giữ lại mồi nhận diện được tối thiểu X% chủng mục tiêu.",
        "pcr_product_size": "Kích thước Sp PCR tối thiểu (bp)",
        "max_bp": "Tối đa (bp)",
        "primer_tm": "Nhiệt độ Tm Mồi (°C)",
        "primer_len": "Độ dài mồi (bp)",
        "ai_integration": "🤖 Tích Hợp AI",
        "enable_ai_toggle": "Kích hoạt Chuyên gia AI (AI Expert) đánh giá mồi",
        "enable_ai_help": "AI sẽ tự động kiểm duyệt chéo các bộ mồi và lọc ra cặp hoàn hảo nhất.",
        "use_blast_toggle": "Tìm kiếm Gene đích (NCBI BLAST)",
        "use_blast_help": "Hệ thống sẽ dịch mã sang Protein để xem vùng bảo tồn thuộc gene độc lực nào.",
        "degenerate_toggle": "🔮 Mồi Thoái Hóa (Bậc cao)",
        "degenerate_help": "Sinh mồi thoái hoá (IUPAC) để bắt được các chủng bị đột biến điểm.",
        "expert_mode": "⚙️ Expert Mode & Hardware (Nâng Cao)",
        "ai_server_config": "**1. Cấu hình AI Server**",
        "llm_model": "LLM Model",
        "bio_params": "**2. Thông số Sinh học Chuyên sâu**",
        "min_conservation": "Bảo tồn K-mer (0-1)",
        "min_copy": "Min Intra-Strain Copy",
        "max_cross": "Max Cross-Reactivity (%)",
        "max_mismatches": "Max Mismatches",
        "perf_algo": "**3. Hiệu năng & Thuật toán**",
        "kmer_stride": "K-mer Stride (Step)",
        "kmer_stride_help": "Tăng step để chạy nhanh hơn nhưng dễ lọt mồi.",
        "max_candidates": "Giới hạn Candidate",
        "cpu_cores": "Số Core CPU (0=Auto)",
        "random_sampling": "**Lấy mẫu ngẫu nhiên (0 = Toàn bộ)**",
        "design_target": "Design Target",
        "design_bg": "Design Background",
        "validate_target": "Validate Target",
        "validate_bg": "Validate Background",
        "settings_saved": "Các cài đặt tự động được lưu.",
        "chat_suggest": "💡 **Gợi ý:** Đây là cách dễ nhất để thiết kế mồi. Chỉ cần nói với AI loài bạn muốn thiết kế, hệ thống sẽ tự động cấu hình toàn bộ các thông số kỹ thuật bên trong.",
        "results_loaded": "🧬 Trạng thái Chuyên gia AI (AI Expert): Kết Quả Đã Nạp!",
        "validation_loaded": "🧬 Trạng thái Chuyên gia AI (AI Expert): Kết Quả Kiểm Chứng Đã Nạp!",
        "results_desc": "Phát hiện kết quả thiết kế mồi thành công tại thư mục <code>{dir}</code>. AI đã tự động nạp dữ liệu mồi/probe và sẵn sàng phân tích chuyên sâu!",
        "validation_desc": "Phát hiện kết quả kiểm chứng mồi (Validation) thành công tại thư mục <code>{dir}</code>. AI đã tự động nạp dữ liệu phân tích PCR In-Silico nâng cao và sẵn sàng tư vấn chuyên sâu!",
        "btn_reset": "🗑️ Reset / Thiết Kế Lại Từ Đầu",
        "error_reset": "Không thể reset file kết quả: ",
        "ai_setup_config": "💡 **AI đã thiết lập xong cấu hình chạy:**",
        "target_inclusion": "🎯 **Mục tiêu (Inclusion)**",
        "bg_exclusion": "🛡️ **Bộ nền loại trừ (Exclusion)**",
        "primers_to_validate_label": "🔍 **Các cặp mồi cần kiểm chứng (Primers to Validate):**",
        "confirm_run": "🚀 XÁC NHẬN & CHẠY PIPELINE",
        "error_email": "Vui lòng điền NCBI Email ở Sidebar bên trái trước khi chạy!",
        "chat_history_finish": "Tôi đã thiết lập xong cấu hình chạy phía trên. Bạn có muốn điều chỉnh thêm bớt loài nào không? Nếu ổn thì bạn bấm nút Xác nhận hoặc gõ 'chạy pipeline' nhé!",
        "chat_finished_design": "Tôi đã đánh giá xong kết quả Pipeline!\n\n**🏆 Assay tốt nhất:** {best}\n**📊 Đánh giá tổng quan:** {verdict}\n\nBạn có muốn đi sâu vào phân tích độ nhạy, đặc hiệu hay tính thoái hóa của bộ mồi này không?",
        "chat_input_placeholder": "Yêu cầu thiết kế mồi (VD: Thiết kế mồi cho S. pneumoniae)",
        "tab_chat_title": "💬 Chuyên gia AI (AI Expert)",
        "tab_auto_title": "🤖 Chế Độ Tự Động (NCBI)",
        "tab_local_title": "📂 Dữ Liệu Local",
        "tab_history_title": "📊 Lịch Sử Chạy & Dashboard",
        "auto_header": "Tải Genomes từ NCBI & Thiết Kế",
        "local_header": "Sử dụng file FASTA trên máy",
        "target_query_label": "Target Query",
        "bg_query_label": "Background Query",
        "min_mb": "Min Mb",
        "max_count": "Max Count",
        "add_target": "➕ Thêm loài Target",
        "add_bg": "➕ Thêm loài Background",
        "btn_start_auto": "🚀 Bắt Đầu Thiết Kế (NCBI)",
        "btn_start_local": "🚀 Bắt Đầu Thiết Kế (Local)",
        "project_name": "Tên Project",
        "local_target_path": "Đường dẫn thư mục Target",
        "local_bg_path": "Đường dẫn thư mục Background",
        "local_out_path": "Thư mục đầu ra (Tùy chọn - Tự động định tuyến)",
        "browse": "📂 Chọn",
        "invalid_target": "❌ Đường dẫn Target không hợp lệ.",
        "history_header": "📊 Lịch Sử Thiết Kế & Analytics Dashboard",
        "history_subheader": "Xem lại, phân tích và so sánh các bộ mồi đã thiết kế trong cơ sở dữ liệu hệ thống.",
        "total_targets": "Tổng số mục tiêu đã thiết kế",
        "total_runs": "Tổng số lượt chạy thành công",
        "latest_run": "Lượt chạy mới nhất",
        "no_history": "🗂️ Chưa có lượt chạy nào được ghi nhận trong thư mục `runs/`.",
        "df_species": "Loài / Gene",
        "df_time": "Thời gian chạy",
        "df_best": "Cặp mồi tốt nhất",
        "df_sensitivity": "Độ nhạy (Sensitivity)",
        "df_specificity": "Độ đặc hiệu (Specificity)",
        "df_target_gene": "Gen mục tiêu (BLAST)",
        "df_candidates": "Số ứng viên",
        "df_folder": "Thư mục",
        "detail_explorer": "### 🔎 Detailed Run Explorer",
        "select_run": "Chọn lượt chạy để phân tích sâu:",
        "detail_run_for": "#### 📂 Detailed results for: **{target}**",
        "detail_path": "📁 *Đường dẫn lưu kết quả:* `{path}`",
        "ai_recommendation": "**🤖 Khuyến nghị Lâm sàng từ AI Expert:** ",
        "best_primer_pair": "🏆 **Cặp mồi đề xuất:** ",
        "overall_evaluation": "📊 **Đánh giá tổng quan:** ",
        "quick_assays_title": "🏆 Xem nhanh Top 5 Cặp mồi thiết kế tốt nhất",
        "quick_validation_title": "📊 Xem nhanh kết quả chẩn đoán PCR In-Silico chi tiết",
        "pcr_expert_header": "🤖 Phân tích Sinh học Phân tử từ AI Expert",
        "pcr_analysis_tm": "**🌡️ Phân tích Tm & GC:**\n",
        "pcr_analysis_balance": "**🎯 Cân bằng Đặc hiệu/Độ nhạy:**\n",
        "pcr_analysis_rec": "**💡 Khuyến nghị tối ưu Lâm sàng:** ",
        "structural_warning": "⚠️ Rủi ro cấu trúc: ",
        "no_assays_found": "Không tìm thấy dữ liệu assay chi tiết.",
        "ncbi_email_disabled": "NCBI Email (Sử dụng ở Sidebar)",
        "sidebar_placeholder": "Dùng mục Sidebar",
        "target_group_header": "#### 🎯 Target Group (Inclusion)",
        "bg_group_header": "#### 🛡️ Background Group (Exclusion)",
        "target_query_num": "Từ khóa Target #{num}",
        "bg_query_num": "Từ khóa Background #{num}",
        "min_mb_num": "Min Mb #{num}",
        "max_count_num": "Max Count #{num}",
        "error_email_valid": "Vui lòng nhập NCBI Email hợp lệ ở Sidebar.",
        "runs_placeholder": "runs/ [Tự động định tuyến]",
        "exec_progress": "📝 Tiến trình thực thi",
        "show_logs": "Hiển thị log chi tiết",
        "pipeline_started": "🚀 Bắt đầu chạy pipeline...",
        "stage1_building": "📦 Khởi tạo cơ sở dữ liệu (Stage 1)...",
        "auto_opt_cycle": "🔄 Chu kỳ tối ưu hóa tự động...",
        "stage2_designing": "🧬 Thiết kế mồi (Stage 2)...",
        "stage3_validation": "💻 Kiểm nghiệm In-Silico (Stage 3)...",
        "stage4_probe": "🔍 Thiết kế Probe (Stage 4)...",
        "stage5_blast": "🌐 Định danh gen BLAST (Stage 5)...",
        "stage6_ai": "🤖 Tham vấn AI Core (Stage 6)...",
        "pipeline_success": "✨ Pipeline hoàn thành xuất sắc!",
        "pipeline_failed": "❌ Chạy Pipeline thất bại. Vui lòng kiểm tra log phía trên.",
        "multiplex_header": "Thiết Kế Multiplex PCR (Đa Tác Nhân)",
        "add_multiplex_target": "➕ Thêm loài Target Multiplex",
        "quick_multiplex_assays_title": "🏆 Xem nhanh Top Multiplex Kits Thiết Kế Tốt Nhất",
        "physical_risk_auditing": "🚨 Đánh Giá Rủi Ro Tương Tác Vật Lý & Kiểm Duyệt Dimer",
        "clinical_advice_multiplex": "🤖 Khuyến Nghị Lâm Sàng & qPCR Wet-lab Của AI Expert",
        "recommended_kit": "🏆 Bộ Kit Đề Xuất",
        "clinical_verdict": "📊 Phê Duyệt Lâm Sàng",
        "tm_balance_gc": "Cân Bằng Tm & GC",
        "spec_analysis": "Phân Tích Độ Đặc Hiệu",
        "wet_lab_guide": "Hướng dẫn Tối ưu Hóa Wet-Lab"
    },
    "en": {
        "title": "🤖 V-Extreme AI Expert",
        "subtitle": "AI-Powered Rational Primer Design System",
        "welcome": (
            "👋 **Hello! I am V-Extreme AI Expert.**\n\n"
            "I am here to help you design PCR primers as easily as possible. You can ask me to:\n\n"
            "1-Design new primers using online NCBI data\n\n"
            "2-Design new primers using offline local data from your directory.\n\n"
            "3-Verify the sensitivity and specificity of a known/custom primer set\n\n"
            "4-Analyze in-silico PCR diagnostics for prior design results\n\n"
            "5-Provide wet-lab troubleshooting and qPCR protocol optimization for designed primers/probes."
        ),
        "lang_selector": "🌐 Select Language",
        "reset_search": "🗑️ Reset Search Fields",
        "essential_settings": "🔬 Essential Settings",
        "ncbi_email": "NCBI Email (Required to download Genomes)",
        "save_email": "💾 Save email for next time",
        "min_sensitivity": "Min Sensitivity (%)",
        "min_sensitivity_help": "Filter out primers that do not cover at least X% of target strains.",
        "pcr_product_size": "Min PCR Product Size (bp)",
        "max_bp": "Max (bp)",
        "primer_tm": "Primer Tm Range (°C)",
        "primer_len": "Primer Length (bp)",
        "ai_integration": "🤖 AI Integration",
        "enable_ai_toggle": "Enable AI Expert Evaluation",
        "enable_ai_help": "The AI will perform strict structural and clinical audits of all designed candidates.",
        "use_blast_toggle": "NCBI BLAST Gene Annotation",
        "use_blast_help": "Translates conserved amplicon to protein to find virulence/structural target genes.",
        "degenerate_toggle": "🔮 Degenerate Primers (Advanced)",
        "degenerate_help": "Allow degenerate IUPAC bases to tolerate single nucleotide polymorphisms (SNPs).",
        "expert_mode": "⚙️ Expert Mode & Hardware Settings",
        "ai_server_config": "**1. AI Server Configuration**",
        "llm_model": "LLM Model",
        "bio_params": "**2. Detailed Biochemical Parameters**",
        "min_conservation": "K-mer Conservation (0-1)",
        "min_copy": "Min Intra-Strain Copy",
        "max_cross": "Max Cross-Reactivity (%)",
        "max_mismatches": "Max Mismatches",
        "perf_algo": "**3. Performance & Stride**",
        "kmer_stride": "K-mer Stride (Step)",
        "kmer_stride_help": "Increase stride for higher design speed, lower stride for absolute sensitivity.",
        "max_candidates": "Max Candidate Candidates",
        "cpu_cores": "CPU Cores (0=Auto)",
        "random_sampling": "**Random Subsampling (0 = Full Database)**",
        "design_target": "Design Target",
        "design_bg": "Design Background",
        "validate_target": "Validate Target",
        "validate_bg": "Validate Background",
        "settings_saved": "Settings saved automatically.",
        "chat_suggest": "💡 **Tip:** This is the easiest way to design primers. Just tell the AI what species you want to target, and it will configure all biophysical parameters automatically.",
        "results_loaded": "🧬 AI Expert Status: Results Loaded!",
        "validation_loaded": "🧬 AI Expert Status: Validation Results Loaded!",
        "results_desc": "Successfully loaded designed primers from <code>{dir}</code>. AI is ready to provide clinical or wet-lab troubleshooting advice!",
        "validation_desc": "Successfully loaded in-silico validation results from <code>{dir}</code>. AI is ready to analyze diagnostic coverage and cross-reactivity!",
        "btn_reset": "🗑️ Reset / Redesign From Scratch",
        "error_reset": "Unable to reset result files: ",
        "ai_setup_config": "💡 **AI has configured the following pipeline settings:**",
        "target_inclusion": "🎯 **Target (Inclusion)**",
        "bg_exclusion": "🛡️ **Background (Exclusion)**",
        "primers_to_validate_label": "🔍 **Primers to Validate:**",
        "confirm_run": "🚀 CONFIRM & RUN PIPELINE",
        "error_email": "Please enter your NCBI Email in the Sidebar first!",
        "chat_history_finish": "I have configured the pipeline parameters. Would you like to adjust or add any species? If it looks good, click Confirm above or type 'run pipeline'!",
        "chat_finished_design": "I have completed the pipeline evaluation!\n\n**🏆 Best Assay:** {best}\n**📊 Overall Verdict:** {verdict}\n\nWould you like to analyze structural warning, Tm differences, or wet-lab cycle conditions for this set?",
        "chat_input_placeholder": "Ask to design primers (e.g. Design primers for S. pneumoniae)",
        "tab_chat_title": "💬 AI Expert",
        "tab_auto_title": "🤖 Automatic Mode (NCBI)",
        "tab_local_title": "📂 Local File Mode",
        "tab_history_title": "📊 History & Dashboard",
        "auto_header": "Download Genomes from NCBI & Design",
        "local_header": "Use Local FASTA Files",
        "target_query_label": "Target Query",
        "bg_query_label": "Background Query",
        "min_mb": "Min Mb",
        "max_count": "Max Count",
        "add_target": "➕ Add Target",
        "add_bg": "➕ Add Background",
        "btn_start_auto": "🚀 Start NCBI Pipeline",
        "btn_start_local": "🚀 Start Local Pipeline",
        "project_name": "Project Name",
        "local_target_path": "Path to Target Folder",
        "local_bg_path": "Path to Background Folder",
        "local_out_path": "Output Folder (Optional)",
        "browse": "📂 Browse",
        "invalid_target": "❌ Invalid Target path.",
        "history_header": "📊 Design History & Analytics Dashboard",
        "history_subheader": "Review, analyze, and cross-compare designed assay kits across the database.",
        "total_targets": "Total Target Species Designed",
        "total_runs": "Total Successful Runs",
        "latest_run": "Latest Target Designed",
        "no_history": "🗂️ No runs found inside the `runs/` database.",
        "df_species": "Species / Gene",
        "df_time": "Run Time",
        "df_best": "Best Assay Set",
        "df_sensitivity": "Sensitivity (%)",
        "df_specificity": "Specificity (%)",
        "df_target_gene": "Target Gene (BLAST)",
        "df_candidates": "Candidates Count",
        "df_folder": "Folder",
        "detail_explorer": "### 🔎 Detailed Run Explorer",
        "select_run": "Select a run to perform a deep-dive analysis:",
        "detail_run_for": "#### 📂 Detailed results for: **{target}**",
        "detail_path": "📁 *Output Directory:* `{path}`",
        "ai_recommendation": "**🤖 Clinical Advice from AI Expert:** ",
        "best_primer_pair": "🏆 **Recommended Assay Kit:** ",
        "overall_evaluation": "📊 **Overall Verdict:** ",
        "quick_assays_title": "🏆 Top 5 Designed Assays Summary",
        "quick_validation_title": "📊 Detailed In-Silico PCR Diagnostic Table",
        "pcr_expert_header": "🤖 Molecular Biophysics Audit by AI Expert",
        "pcr_analysis_tm": "**🌡️ Tm & GC Composition:**\n",
        "pcr_analysis_balance": "**🎯 Inclusivity & Exclusivity Profile:**\n",
        "pcr_analysis_rec": "**💡 Clinical Optimization & Troubleshooting:** ",
        "structural_warning": "⚠️ Physical warning: ",
        "no_assays_found": "No detailed assays found.",
        "ncbi_email_disabled": "NCBI Email (Use Sidebar)",
        "sidebar_placeholder": "Use Sidebar",
        "target_group_header": "#### 🎯 Target Group (Inclusion)",
        "bg_group_header": "#### 🛡️ Background Group (Exclusion)",
        "target_query_num": "Target Query #{num}",
        "bg_query_num": "Background Query #{num}",
        "min_mb_num": "Min Mb #{num}",
        "max_count_num": "Max Count #{num}",
        "error_email_valid": "Please enter a valid NCBI email in the Sidebar.",
        "runs_placeholder": "runs/ [Auto Routed]",
        "exec_progress": "📝 Execution Progress",
        "show_logs": "Show detailed logs",
        "pipeline_started": "🚀 Pipeline started...",
        "stage1_building": "📦 Building Datasets (Stage 1)...",
        "auto_opt_cycle": "🔄 Auto-Optimization Cycle...",
        "stage2_designing": "🧬 Designing Primers (Stage 2)...",
        "stage3_validation": "💻 In-Silico Validation (Stage 3)...",
        "stage4_probe": "🔍 Probe Design (Stage 4)...",
        "stage5_blast": "🌐 BLAST Annotation (Stage 5)...",
        "stage6_ai": "🤖 AI Core Consultation (Stage 6)...",
        "pipeline_success": "✨ Pipeline Completed Successfully!",
        "pipeline_failed": "❌ Pipeline Failed. Check logs above.",
        "multiplex_header": "Multiplex PCR Kit Design",
        "add_multiplex_target": "➕ Add Multiplex Target Species",
        "quick_multiplex_assays_title": "🏆 Top Multiplex Kits Summary",
        "physical_risk_auditing": "🚨 Physical Biophysical Interactions & Dimer Auditing",
        "clinical_advice_multiplex": "🤖 AI Expert Clinical Advice & qPCR Wet-lab Recommendations",
        "recommended_kit": "RECOMMENDED KIT",
        "clinical_verdict": "CLINICAL VERDICT",
        "tm_balance_gc": "Tm Balance & GC Composition",
        "spec_analysis": "Specificity & Sensitivity Balance",
        "wet_lab_guide": "Wet-lab Protocol Suggestions"
    }
}

def tr(key):
    lang = st.session_state.get("language", "vi")
    return T[lang].get(key, key)

st.set_page_config(page_title="Rational Primer Design", page_icon="🧬", layout="wide")

# --- TKINTER SAFETY CHECK ---
SYSTEM_OS = platform.system()
HAS_TK = False

# Spawning a quick subprocess to check if tkinter is available on this system,
# avoiding any initialization of Tk/Tcl in Streamlit's worker threads.
try:
    import subprocess
    import sys
    res = subprocess.run([sys.executable, "-c", "import tkinter"], capture_output=True, timeout=2)
    if res.returncode == 0:
        HAS_TK = True
except Exception:
    HAS_TK = False

# --- 2. INITIALIZE SESSION STATE ---
if "language" not in st.session_state:
    st.session_state["language"] = "vi"

lang = st.session_state["language"]

if "reset_id" not in st.session_state: st.session_state["reset_id"] = 0 

if "email_val" not in st.session_state:
    saved_email = ""
    st.session_state["save_email_cb"] = False
    if os.path.exists("config/user_settings.json"):
        try:
            with open("config/user_settings.json", "r") as f:
                saved_email = json.load(f).get("saved_email", "")
                st.session_state["save_email_cb"] = bool(saved_email)
        except Exception:
            pass
    st.session_state["email_val"] = saved_email
if "path_t_val" not in st.session_state: st.session_state["path_t_val"] = ""
if "path_b_val" not in st.session_state: st.session_state["path_b_val"] = ""
if "l_out_val" not in st.session_state: st.session_state["l_out_val"] = os.path.join(os.getcwd(), "results_local")
if "out_dir_val" not in st.session_state: st.session_state["out_dir_val"] = "results_auto"
if "auto_proj_val" not in st.session_state: st.session_state["auto_proj_val"] = "Auto_Run_01"
if "local_proj_val" not in st.session_state: st.session_state["local_proj_val"] = "Local_Run_01"

if "chat_history" not in st.session_state:
    st.session_state["chat_history"] = [{"role": "assistant", "content": T[lang]["welcome"]}]
else:
    # Update welcome message dynamically if it's the only message and language changed
    if len(st.session_state["chat_history"]) == 1 and st.session_state["chat_history"][0]["role"] == "assistant":
        if st.session_state["chat_history"][0]["content"] in [T["vi"]["welcome"], T["en"]["welcome"]]:
            st.session_state["chat_history"][0]["content"] = T[lang]["welcome"]

if "copilot_json" not in st.session_state:
    st.session_state["copilot_json"] = None

if "target_list" not in st.session_state:
    st.session_state["target_list"] = [{"query": "", "size": 0.0, "count": 500}]
if "bg_list" not in st.session_state:
    st.session_state["bg_list"] = [{"query": "", "size": 0.0, "count": 50}]
if "multiplex_target_list" not in st.session_state:
    st.session_state["multiplex_target_list"] = [{"query": "", "size": 0.0, "count": 500}, {"query": "", "size": 0.0, "count": 500}]
if "is_multiplex" not in st.session_state:
    st.session_state["is_multiplex"] = False
if "multiplex_type" not in st.session_state:
    st.session_state["multiplex_type"] = "qPCR"
if "is_local_multiplex" not in st.session_state:
    st.session_state["is_local_multiplex"] = False
if "local_target_list" not in st.session_state:
    st.session_state["local_target_list"] = ["", ""]

# --- 3. HELPER FUNCTIONS ---

def generate_run_output_dir(target_name):
    """Tạo đường dẫn thư mục lưu kết quả chuẩn hóa theo tên loài/gene và thời gian chạy."""
    import re
    from datetime import datetime
    
    cleaned = re.sub(r'[^a-zA-Z0-9_\s-]', '', target_name)
    cleaned = cleaned.strip().replace(' ', '_').lower()
    if not cleaned:
        cleaned = "unknown_target"
        
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    runs_dir = os.path.join(os.getcwd(), "runs")
    os.makedirs(runs_dir, exist_ok=True)
    
    target_dir = os.path.join(runs_dir, cleaned)
    os.makedirs(target_dir, exist_ok=True)
    
    run_dir = os.path.join(target_dir, f"run_{timestamp}")
    os.makedirs(run_dir, exist_ok=True)
    return run_dir

def reset_app():
    st.session_state["path_t_val"] = ""
    st.session_state["path_b_val"] = ""
    st.session_state["l_out_val"] = os.path.join(os.getcwd(), "results_local")
    st.session_state["out_dir_val"] = "results_auto"
    st.session_state["auto_proj_val"] = "Auto_Run_01"
    st.session_state["local_proj_val"] = "Local_Run_01"
    st.session_state["target_list"] = [{"query": "", "size": 0.0, "count": 500}]
    st.session_state["bg_list"] = [{"query": "", "size": 0.0, "count": 50}]
    st.session_state["multiplex_target_list"] = [{"query": "", "size": 0.0, "count": 500}, {"query": "", "size": 0.0, "count": 500}]
    st.session_state["is_multiplex"] = False
    st.session_state["multiplex_type"] = "qPCR"
    st.session_state["is_local_multiplex"] = False
    st.session_state["local_target_list"] = ["", ""]
    st.session_state["reset_id"] += 1

def select_folder_callback(session_key):
    if not HAS_TK:
        return
    try:
        import subprocess
        import sys
        
        # Spawning a separate subprocess to display the Tkinter folder browse dialog.
        # This guarantees it runs on the main thread of the subprocess, preventing
        # macOS "notifier not initialized" / Abort trap 6 crashes inside Streamlit's worker threads.
        code = (
            "import tkinter as tk\n"
            "from tkinter import filedialog\n"
            "root = tk.Tk()\n"
            "root.withdraw()\n"
            "root.wm_attributes('-topmost', 1)\n"
            "path = filedialog.askdirectory()\n"
            "root.destroy()\n"
            "if path:\n"
            "    print(path)"
        )
        res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
        folder_path = res.stdout.strip()
        if folder_path:
            st.session_state[session_key] = folder_path
    except Exception as e:
        st.error(f"Error opening folder dialog: {e}")

# --- 4. CSS STYLING ---
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&display=swap');
    
    /* Global Styling */
    html, body, [class*="css"], .stApp {
        font-family: 'Outfit', sans-serif !important;
    }
    
    .reportview-container { margin-top: -2em; }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stDeployButton {display:none;}
    .stButton button { margin-top: 0px; }
    
    div.stCode > div > pre {
        max-height: 500px;
        overflow-y: auto !important;
        white-space: pre-wrap !important;
    }
    
    /* Modern Glassmorphic Tabs */
    div[data-testid="stTabBar"] {
        background-color: #f8fafc;
        padding: 8px;
        border-radius: 12px;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.01);
        border: 1px solid #e2e8f0;
        gap: 4px;
        margin-bottom: 20px;
    }
    button[data-testid="stMarkdownContainer"] p {
        font-size: 15px !important;
        font-weight: 600 !important;
        margin: 0 !important;
    }
    div[data-testid="stTabBar"] button {
        border-radius: 8px !important;
        padding: 8px 18px !important;
        transition: all 0.2s ease-in-out !important;
        border: none !important;
        background: transparent !important;
    }
    div[data-testid="stTabBar"] button[aria-selected="true"] {
        background-color: #ffffff !important;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05), 0 2px 4px -1px rgba(0, 0, 0, 0.03) !important;
        border-bottom: 2px solid #3b82f6 !important;
        color: #1e3a8a !important;
    }
    
    /* Metric Widgets Styling */
    div[data-testid="stMetricValue"] {
        font-size: 26px !important;
        font-weight: 700 !important;
        color: #10b981 !important;
    }
    div[data-testid="stMetricLabel"] {
        font-weight: 600 !important;
        color: #10b981 !important;
    }
    
    /* Modernized Expanders */
    .streamlit-expanderHeader {
        background-color: #ffffff !important;
        border-radius: 8px !important;
        border: 1px solid #e2e8f0 !important;
        box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.01) !important;
        transition: border-color 0.2s, box-shadow 0.2s !important;
    }
    .streamlit-expanderHeader:hover {
        border-color: #3b82f6 !important;
        box-shadow: 0 4px 6px -1px rgba(0,0,0,0.02) !important;
    }
    
    /* Premium Title Header Container */
    .custom-header {
        text-align: center;
        padding: 25px 20px;
        background: linear-gradient(135deg, #1e3a8a 0%, #0f172a 100%);
        border-radius: 16px;
        box-shadow: 0 10px 25px -5px rgba(15, 23, 42, 0.2), 0 8px 10px -6px rgba(15, 23, 42, 0.1);
        border: 1px solid rgba(255, 255, 255, 0.08);
        color: white;
        margin-bottom: 25px;
    }
    .custom-title {
        font-family: 'Outfit', sans-serif !important;
        font-weight: 700 !important;
        font-size: 2.2rem !important;
        letter-spacing: -0.04rem !important;
        background: linear-gradient(90deg, #38bdf8 0%, #3b82f6 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin: 0 !important;
    }
    .custom-subtitle {
        font-family: 'Outfit', sans-serif !important;
        font-weight: 500 !important;
        font-size: 1.15rem !important;
        color: #e2e8f0 !important;
        margin-top: 6px !important;
        margin-bottom: 2px !important;
    }
</style>
""", unsafe_allow_html=True)

st.markdown(f"""
<div class="custom-header">
    <h1 class="custom-title">{tr("title")}</h1>
    <p class="custom-subtitle">{tr("subtitle")}</p>
</div>
""", unsafe_allow_html=True)
st.markdown("<br>", unsafe_allow_html=True)

# --- 5. SIDEBAR ---
with st.sidebar:
    # 🌐 CHỌN NGÔN NGỮ / SELECT LANGUAGE
    lang_choice = st.radio("🌐 Chọn ngôn ngữ / Language", ["Tiếng Việt", "English"], index=0 if st.session_state.get("language", "vi") == "vi" else 1, horizontal=True)
    lang = "vi" if lang_choice == "Tiếng Việt" else "en"
    if st.session_state.get("language") != lang:
        st.session_state["language"] = lang
        st.rerun()
        
    if st.button(tr("reset_search"), on_click=reset_app):
        st.rerun()
    
    st.divider()
    
    st.header(tr("essential_settings"))
    email_val = st.text_input(tr("ncbi_email"), value=st.session_state["email_val"], placeholder="email@example.com")
    st.session_state["email_val"] = email_val
    save_email = st.checkbox(tr("save_email"), value=st.session_state.get("save_email_cb", False))
    st.session_state["save_email_cb"] = save_email
    
    if save_email and email_val:
        os.makedirs("config", exist_ok=True)
        with open("config/user_settings.json", "w") as f:
            json.dump({"saved_email": email_val}, f)
    elif not save_email and os.path.exists("config/user_settings.json"):
        try: os.remove("config/user_settings.json")
        except: pass

    st.subheader(tr("primer_tm"))
    min_sens = st.slider(tr("min_sensitivity"), 50.0, 100.0, 100.0, 0.1, help=tr("min_sensitivity_help"))
    
    c1, c2 = st.columns(2)
    prod_min = c1.number_input(tr("pcr_product_size"), value=90)
    prod_max = c2.number_input(tr("max_bp"), value=200)
    
    primer_tm_range = st.slider(tr("primer_tm"), 40.0, 75.0, (55.0, 65.0), 0.5)
    tm_min, tm_max = primer_tm_range
    p_len_range = st.slider(tr("primer_len"), 15, 30, (18, 22))

    st.divider()
    st.header(tr("ai_integration"))
    enable_ai = st.toggle(tr("enable_ai_toggle"), value=True, help=tr("enable_ai_help"), key="enable_ai_cb")
    use_blast = st.toggle(tr("use_blast_toggle"), value=True, help=tr("use_blast_help"), key="use_blast_cb")
    deg_primers = st.toggle(tr("degenerate_toggle"), value=False, help=tr("degenerate_help"), key="degenerate_cb")

    with st.expander(tr("expert_mode"), expanded=False):
        st.markdown(tr("ai_server_config"))
        ai_base_url = st.text_input("Local API Base URL", value="http://localhost:11434/v1", key="ai_base_url_val")
        @st.cache_data(ttl=10, show_spinner=False)
        def fetch_models(base_url):
            try:
                from openai import OpenAI
                client = OpenAI(api_key="sk-local", base_url=base_url)
                return [m.id for m in client.models.list().data]
            except: return []
        available_models = fetch_models(ai_base_url)
        if available_models:
            ai_model = st.selectbox(tr("llm_model"), available_models, key="ai_model_val")
        else:
            ai_model = st.text_input(tr("llm_model"), value="llama3", key="ai_model_val")

        st.markdown("---")
        st.markdown(tr("bio_params"))
        min_cons = st.slider(tr("min_conservation"), 0.5, 1.0, 1.0, 0.01)
        min_copy = st.slider(tr("min_copy"), 1, 10, 2)
        max_xr = st.slider(tr("max_cross"), 0.0, 10.0, 1.0, 0.1)
        max_mm = st.number_input(tr("max_mismatches"), value=2)

        st.markdown("---")
        st.markdown(tr("perf_algo"))
        kmer_step = st.slider(tr("kmer_stride"), 1, 10, 1, help=tr("kmer_stride_help"))
        max_cand = st.slider(tr("max_candidates"), 10, 2000, 10, 10)
        cpu = st.number_input(tr("cpu_cores"), value=0)
        st.markdown(tr("random_sampling"))
        samp_dt = st.number_input(tr("design_target"), value=0)
        samp_db = st.number_input(tr("design_bg"), value=100)
        samp_vt = st.number_input(tr("validate_target"), value=0)
        samp_vb = st.number_input(tr("validate_bg"), value=200)

    st.markdown("---")
    st.info(tr("settings_saved"))

# --- 6. CORE LOGIC ---

def save_params():
    params = {
        "target_mode": "universal",
        "min_sensitivity": min_sens,
        "design_min_conservation": min_cons,
        "min_intra_strain_copy": min_copy,
        "kmer_step": kmer_step, # <--- LƯU THAM SỐ STEP VÀO FILE CONFIG
        "primer_tm_min": tm_min,
        "primer_tm_max": tm_max,
        "primer_opt_tm": (tm_min + tm_max) / 2,
        "validation_max_cross_reactivity": max_xr,
        "product_size_min": prod_min,
        "product_size_max": prod_max,
        "primer_length_min": p_len_range[0],
        "primer_length_max": p_len_range[1],
        "max_mismatch": max_mm,
        "cpu_cores": cpu,
        "design_target_sampling_size": samp_dt,
        "design_background_sampling_size": samp_db,
        "validation_target_sampling_size": samp_vt,
        "validation_background_sampling_size": samp_vb,
        "design_max_candidates": max_cand,
        "enable_blast": use_blast,
        "degenerate_primers": deg_primers
    }
    os.makedirs("config", exist_ok=True)
    path = "config/gui_params.json"
    with open(path, "w") as f:
        json.dump(params, f, indent=4)
    return path

def run_pipeline(cmd, out_dir):
    import shlex

    def command_chain(command):
        if isinstance(command, (list, tuple)):
            return [list(command)]
        tokens = shlex.split(command)
        commands, current = [], []
        for token in tokens:
            if token == "&&":
                if current:
                    commands.append(current)
                    current = []
            else:
                current.append(token)
        if current:
            commands.append(current)
        return commands

    st.divider()
    st.header(tr("exec_progress"))
    
    # Progress Bar Interface
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    with st.expander(tr("show_logs"), expanded=True):
        log_area = st.empty()
    logs = []

    my_env = os.environ.copy()
    my_env["PYTHONIOENCODING"] = "utf-8"
    
    status_text.info(tr("pipeline_started"))

    return_code = 1
    for command_args in command_chain(cmd):
        process = subprocess.Popen(
            command_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=False,
            text=True,
            bufsize=1,
            env=my_env,
            encoding="utf-8"
        )

        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line:
                clean_line = line.strip()
                logs.append(clean_line)

                # Simple Log Parser for Progress
                if "[STAGE 1]" in clean_line:
                    progress_bar.progress(10)
                    status_text.info(tr("stage1_building"))
                elif "CYCLE" in clean_line:
                    progress_bar.progress(20)
                    status_text.info(tr("auto_opt_cycle"))
                elif "[STAGE 2]" in clean_line:
                    progress_bar.progress(30)
                    status_text.info(tr("stage2_designing"))
                elif "[Phase 1]" in clean_line:
                    progress_bar.progress(40)
                elif "[Phase 2]" in clean_line:
                    progress_bar.progress(50)
                elif "[Phase 3]" in clean_line:
                    progress_bar.progress(60)
                elif "[STAGE 3]" in clean_line:
                    progress_bar.progress(70)
                    status_text.info(tr("stage3_validation"))
                elif "[STAGE 4]" in clean_line:
                    progress_bar.progress(85)
                    status_text.info(tr("stage4_probe"))
                elif "[STAGE 5]" in clean_line:
                    progress_bar.progress(95)
                    status_text.info(tr("stage5_blast"))
                elif "[STAGE 6]" in clean_line:
                    progress_bar.progress(98)
                    status_text.info(tr("stage6_ai"))

                # Limit logs to last 150 lines to keep UI snappy
                log_area.code("\n".join(logs[-150:]), language="bash")

        return_code = process.returncode
        if return_code != 0:
            break

    if return_code == 0:
        progress_bar.progress(100)
        status_text.success(tr("pipeline_success"))
        st.balloons()
        
        # --- POST-PROCESSING FOR CUSTOM PRIMER VALIDATION ---
        val_csv = os.path.join(out_dir, "PCR_Advanced_Report.csv")
        if os.path.exists(val_csv) and not os.path.exists(os.path.join(out_dir, "ai_expert_report.json")):
            try:
                import pandas as pd
                df = pd.read_csv(val_csv)
                target_df = df[df["group"].str.lower() == "target"]
                bg_df = df[df["group"].str.lower() == "background"]
                
                total_target = len(target_df["id"].unique()) if not target_df.empty else 0
                pos_target = len(target_df[target_df["status"].str.lower() == "positive"]["id"].unique()) if not target_df.empty else 0
                sens = (pos_target / total_target * 100) if total_target > 0 else 0.0
                
                total_bg = len(bg_df["id"].unique()) if not bg_df.empty else 0
                pos_bg = len(bg_df[bg_df["status"].str.lower() == "positive"]["id"].unique()) if not bg_df.empty else 0
                spec = ((total_bg - pos_bg) / total_bg * 100) if total_bg > 0 else 100.0
                
                markers = df["marker"].unique()
                marker_str = ", ".join(markers) if len(markers) > 0 else "N/A"
                
                if lang == "en":
                    verdict = f"Sensitivity {sens:.1f}% | Specificity {spec:.1f}%"
                    tm_desc = f"Evaluating primer set [{marker_str}] across {total_target} Target strains and {total_bg} Background strains."
                    spec_sens_desc = f"PCR positive on {pos_target}/{total_target} target strains and cross-reactive on {pos_bg}/{total_bg} background strains."
                    
                    if sens >= 95.0 and spec >= 95.0:
                        rec = "The primer set achieves outstanding clinical sensitivity and specificity (>95%). Fully qualified for real wet-lab qPCR deployment."
                    elif sens < 90.0:
                        rec = "Low sensitivity (<90%). Forward/Reverse primers have mismatch alignments on multiple strains. Recommend designing degenerate primers or shifting to another conserved region."
                    elif spec < 90.0:
                        rec = "Poor specificity (<90%). Severe cross-reactivity with background species detected. Recommend altering primer sequences or designing a high-specificity TaqMan probe."
                    else:
                        rec = "The primer set has moderate performance. Further optimize cycle profiles or adjust Mg2+ concentration in the wet-lab to improve efficiency."
                else:
                    verdict = f"Độ nhạy {sens:.1f}% | Độ đặc hiệu {spec:.1f}%"
                    tm_desc = f"Đánh giá bộ mồi [{marker_str}] trên {total_target} chủng Target và {total_bg} chủng Background."
                    spec_sens_desc = f"PCR dương tính trên {pos_target}/{total_target} chủng mục tiêu và phản ứng chéo trên {pos_bg}/{total_bg} chủng ngoại nhóm."
                    
                    if sens >= 95.0 and spec >= 95.0:
                        rec = "Bộ mồi đạt độ nhạy và đặc hiệu xuất sắc lâm sàng (>95%). Hoàn toàn đủ tiêu chuẩn để triển khai thử nghiệm thực tế (wet-lab) qPCR."
                    elif sens < 90.0:
                        rec = "Độ nhạy thấp (<90%). Mồi Forward/Reverse bị lệch mismatch ở nhiều chủng. Khuyến nghị thiết kế lại mồi thoái hóa hoặc đổi vùng bảo tồn khác."
                    elif spec < 90.0:
                        rec = "Độ đặc hiệu kém (<90%). Có phản ứng chéo nghiêm trọng với nhóm ngoại nền (background species). Khuyến nghị thay đổi trình tự mồi hoặc thiết kế probe TaqMan đặc hiệu cao ở giữa."
                    else:
                        rec = "Bộ mồi có hiệu năng ở mức trung bình. Cần tối ưu hóa thêm chu kỳ nhiệt hoặc tăng nồng độ muối Mg2+ ở wet-lab để tăng hiệu suất."
                
                report = {
                    "best_assay_name": marker_str,
                    "overall_verdict": verdict,
                    "tm_analysis": tm_desc,
                    "specificity_sensitivity_balance": spec_sens_desc,
                    "clinical_recommendation": rec,
                    "structural_risks": []
                }
                
                with open(os.path.join(out_dir, "ai_expert_report.json"), "w", encoding="utf-8") as rf:
                    json.dump(report, rf, ensure_ascii=False, indent=4)
            except Exception as e:
                report = {
                    "best_assay_name": "Validation Run",
                    "overall_verdict": "Completed",
                    "tm_analysis": "In-silico PCR execution finished successfully.",
                    "specificity_sensitivity_balance": "Results analyzed.",
                    "clinical_recommendation": f"Lỗi tạo báo cáo tự động: {e}" if lang == "vi" else f"Report generation error: {e}",
                    "structural_risks": []
                }
                with open(os.path.join(out_dir, "ai_expert_report.json"), "w", encoding="utf-8") as rf:
                    json.dump(report, rf, ensure_ascii=False, indent=4)

        report_path = os.path.join(out_dir, "ai_expert_report.json")
        if os.path.exists(report_path):
            with open(report_path, "r", encoding="utf-8") as rf:
                report = json.load(rf)
            
            # Push to Chat History for Interactive AI Evaluator
            if "chat_history" in st.session_state:
                st.session_state["chat_history"].append({
                    "role": "assistant",
                    "content": tr("chat_finished_design").format(best=report.get('best_assay_name', 'N/A'), verdict=report.get('overall_verdict', 'N/A'))
                })
                st.session_state["last_expert_report"] = report
            
            st.divider()
            
            multiplex_csv = os.path.join(out_dir, "MULTIPLEX_KITS.csv")
            if os.path.exists(multiplex_csv):
                render_multiplex_details(out_dir, lang)
            else:
                st.subheader(tr("pcr_expert_header"))
                
                c1, c2 = st.columns([1, 2])
                with c1:
                    st.metric(label=tr("best_primer_pair").replace(":", "").strip(), value=report.get('best_assay_name', 'N/A'))
                    st.metric(label=tr("overall_evaluation").replace(":", "").strip(), value=report.get('overall_verdict', 'N/A'))
                    if 'confidence_score' in report:
                        st.metric(label="💡 " + ("ĐỘ TỰ TIN (AI)" if lang == "vi" else "AI CONFIDENCE"), value=f"{report['confidence_score']}%")
                with c2:
                    st.info(f"{tr('pcr_analysis_tm')}{report.get('tm_analysis', '')}")
                    st.success(f"{tr('pcr_analysis_balance')}{report.get('specificity_sensitivity_balance', '')}")
                    if 'decision_justification' in report:
                        st.warning(f"**Lập luận thuật toán (Algorithm Justification):** {report['decision_justification']}")
                    
                if report.get('structural_risks'):
                    for w in report.get('structural_risks'):
                        st.warning(f"{tr('structural_warning')}{w}")
                
                st.markdown(f"{tr('pcr_analysis_rec')}{report.get('clinical_recommendation', '')}")
    else:
        status_text.error(tr("pipeline_failed"))

def execute_proposal_pipeline(cfg, out_dir_run=None):
    import re
    import json
    import os
    import shutil
    from datetime import datetime
    from rational_design.fetcher import SequenceFetcher
    
    param_file = save_params()
    
    # Override sampling sizes if AI provided them dynamically based on RAM
    if "design_target_sampling_size" in cfg or "design_background_sampling_size" in cfg:
        with open(param_file, "r") as f:
            p = json.load(f)
        if "design_target_sampling_size" in cfg:
            p["design_target_sampling_size"] = cfg["design_target_sampling_size"]
        if "design_background_sampling_size" in cfg:
            p["design_background_sampling_size"] = cfg["design_background_sampling_size"]
        with open(param_file, "w") as f:
            json.dump(p, f, indent=4)
            
    enable_ai = st.session_state.get("enable_ai_cb", True)
    ai_base_url = st.session_state.get("ai_base_url_val", "http://localhost:11434/v1")
    ai_model = st.session_state.get("ai_model_val", "llama3")
    lang = st.session_state.get("language", "vi")
    email = st.session_state.get("email_val", "")
    
    if cfg.get("action") == "propose_multiplex":
        targets = cfg.get("targets", [])
        backgrounds = cfg.get("background", [])
        
        # 1. Determine run folder
        if not out_dir_run:
            target_names = [t["query"] for t in targets]
            proj_name = "multiplex_" + "_".join([re.sub(r'[^a-zA-Z0-9]', '', n.replace(" ", "_").lower()) for n in target_names])
            if len(proj_name) > 60:
                proj_name = proj_name[:60]
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            out_dir_run = os.path.join(os.getcwd(), "runs", proj_name, f"run_{timestamp}")
        
        os.makedirs(out_dir_run, exist_ok=True)
        
        st.info("🔄 Đang chuẩn bị tải dữ liệu multiplex NCBI..." if lang == "vi" else "🔄 Preparing to fetch NCBI multiplex genomes...")
        
        # 2. Download shared background once
        shared_bg_dir = os.path.join(out_dir_run, "0_raw_data", "shared_bg")
        os.makedirs(shared_bg_dir, exist_ok=True)
        fetcher = SequenceFetcher(email)
        
        b_conf = {f"b{bidx+1}": [b["query"], b.get("size", 0.0), b.get("count", 50), b.get("type", "genome")] for bidx, b in enumerate(backgrounds) if b["query"]}
        if b_conf:
            st.info(f"🧬 Đang tải {len(b_conf)} nhóm exclusion background dùng chung..." if lang == "vi" else f"🧬 Fetching {len(b_conf)} shared exclusion background groups...")
            fetcher.fetch_and_save_all(b_conf, shared_bg_dir)
        
        # 3. Download individual target genomes once
        target_dirs = []
        for idx, t in enumerate(targets):
            t_dir = os.path.join(out_dir_run, "0_raw_data", f"t{idx+1}")
            target_dirs.append(t_dir)
            os.makedirs(t_dir, exist_ok=True)
            
            t_conf = {f"t1": [t["query"], t.get("size", 0.0), t.get("count", 500), t.get("type", "genome")]}
            st.info(f"🧬 Đang tải bộ genome target #{idx+1}: {t['query']}..." if lang == "vi" else f"🧬 Fetching target #{idx+1} genome: {t['query']}...")
            fetcher.fetch_and_save_all(t_conf, t_dir)
            
        # 4. Construct cross-pathogen mutual exclusion backgrounds
        # For target i, the background is: shared_bg + targets j (j != i)
        cmd_parts = []
        target_folders = []
        
        for idx, t in enumerate(targets):
            target_folder = os.path.join(out_dir_run, f"t{idx+1}")
            target_folders.append(target_folder)
            os.makedirs(target_folder, exist_ok=True)
            
            # Proactively save target metadata for multiplex analyzer
            try:
                with open(os.path.join(target_folder, "metadata.json"), "w", encoding="utf-8") as mf:
                    json.dump({"target_name": t["query"]}, mf, ensure_ascii=False, indent=4)
            except Exception:
                pass
            
            # Construct a virtual background directory for target i
            virtual_bg_dir = os.path.join(out_dir_run, "0_raw_data", f"bg_for_t{idx+1}")
            os.makedirs(virtual_bg_dir, exist_ok=True)
            
            # Copy shared backgrounds
            if os.path.exists(shared_bg_dir):
                for item in os.listdir(shared_bg_dir):
                    s_path = os.path.join(shared_bg_dir, item)
                    if os.path.isfile(s_path):
                        shutil.copy2(s_path, os.path.join(virtual_bg_dir, item))
                        
            # Copy other target genomes as cross-reactivity backgrounds
            for jdx, other_t_dir in enumerate(target_dirs):
                if jdx != idx:
                    if os.path.exists(other_t_dir):
                        for item in os.listdir(other_t_dir):
                            # Prefix with t{jdx+1}_ to prevent naming collisions
                            s_path = os.path.join(other_t_dir, item)
                            if os.path.isfile(s_path):
                                shutil.copy2(s_path, os.path.join(virtual_bg_dir, f"target_bg_t{jdx+1}_{item}"))
            
            # Write t_conf config file
            t_conf = {f"t1": [t["query"], t.get("size", 0.0), t.get("count", 500), t.get("type", "genome")]}
            t_conf_path = f"config/t_conf_target{idx+1}.json"
            os.makedirs("config", exist_ok=True)
            with open(t_conf_path, "w") as f:
                json.dump(t_conf, f)
                
            # Execute pipeline offline on the pre-fetched local target and background folders!
            cmd_part = f'"{sys.executable}" -m rational_design.cli pipeline --out "{target_folder}" --local_target "{target_dirs[idx]}" --local_bg "{virtual_bg_dir}" --params "{param_file}" --language "{lang}"'
            if enable_ai:
                cmd_part += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
            cmd_parts.append(cmd_part)
            
        assay_type = cfg.get("assay_type") or st.session_state.get("multiplex_type", "qPCR")
        folders_str = " ".join([f'"{f}"' for f in target_folders])
        
        # shared_primer_context.json enables sequential cross-target gating:
        # T1 creates it on success, T2 reads T1's accepted primers and runs Gate 4 before accepting, etc.
        shared_context_file = os.path.join(out_dir_run, "shared_primer_context.json")
        
        # Inject --shared_context into every per-target pipeline command
        cmd_parts_with_ctx = []
        for cp in cmd_parts:
            if "--local_target" in cp:  # Only per-target pipeline commands
                cp += f' --shared_context "{shared_context_file}"'
            cmd_parts_with_ctx.append(cp)
        
        cmd_multi = f'"{sys.executable}" -m rational_design.cli multiplex_analyze --folders {folders_str} --out "{out_dir_run}" --language "{lang}" --assay_type "{assay_type}"'
        if enable_ai:
            cmd_multi += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
        cmd_parts_with_ctx.append(cmd_multi)
        
        cmd = " && ".join(cmd_parts_with_ctx)

        
    elif cfg.get("action") == "propose_validation":
        import csv
        os.makedirs("config", exist_ok=True)
        csv_path = "config/primers_to_validate.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["name", "fwd", "rev"])
            for p in cfg.get("primers", []):
                writer.writerow([p["name"], p["fwd"], p["rev"]])
                
        t_conf = {f"t{idx+1}": [it["query"], it["size"], it["count"], it.get("type", "genome")] for idx, it in enumerate(cfg.get("target", []))}
        b_conf = {f"b{idx+1}": [it["query"], it["size"], it["count"], it.get("type", "genome")] for idx, it in enumerate(cfg.get("background", []))}
        with open("config/t_conf.json", "w") as f: json.dump(t_conf, f)
        with open("config/b_conf.json", "w") as f: json.dump(b_conf, f)
        
        target_species_name = cfg.get("target", [])[0]["query"] if cfg.get("target", []) else "validation_run"
        out_dir_run = generate_run_output_dir(target_species_name)
        
        fetch_cmd = f'"{sys.executable}" -c "import json; from rational_design.fetcher import SequenceFetcher; f=SequenceFetcher(\'{st.session_state["email_val"]}\'); f.fetch_and_save_all(json.load(open(\'config/t_conf.json\')), \'{out_dir_run}/0_raw_data/target\'); f.fetch_and_save_all(json.load(open(\'config/b_conf.json\')), \'{out_dir_run}/0_raw_data/background\')"'
        pcr_cmd = f'"{sys.executable}" -m rational_design.cli validate_primers -c config/primers_to_validate.csv -t "{out_dir_run}/0_raw_data/target" -b "{out_dir_run}/0_raw_data/background" -o "{out_dir_run}/PCR_Advanced_Report.csv"'
        cmd = f"{fetch_cmd} && {pcr_cmd}"
    elif cfg.get("action") == "propose_local_design":
        l_target = cfg.get("local_target", "")
        l_bg = cfg.get("local_background", "")
        target_species_name = os.path.basename(l_target)
        out_dir_run = generate_run_output_dir(target_species_name)
        cmd = f'"{sys.executable}" -m rational_design.cli pipeline --out "{out_dir_run}" --local_target "{l_target}" --local_bg "{l_bg}" --params "{param_file}" --language "{lang}"'
        if enable_ai:
            cmd += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
    else:
        t_conf = {f"t{idx+1}": [it["query"], it["size"], it["count"], it.get("type", "genome")] for idx, it in enumerate(cfg.get("target", []))}
        b_conf = {f"b{idx+1}": [it["query"], it["size"], it["count"], it.get("type", "genome")] for idx, it in enumerate(cfg.get("background", []))}
        with open("config/t_conf.json", "w") as f: json.dump(t_conf, f)
        with open("config/b_conf.json", "w") as f: json.dump(b_conf, f)
        
        target_species_name = cfg.get("target", [])[0]["query"] if cfg.get("target", []) else "unknown"
        out_dir_run = generate_run_output_dir(target_species_name)
        
        cmd = f'"{sys.executable}" -m rational_design.cli pipeline --out "{out_dir_run}" --email "{st.session_state["email_val"]}" --target_config "config/t_conf.json" --bg_config "config/b_conf.json" --params "{param_file}" --language "{lang}"'
        if enable_ai:
            cmd += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
            
    run_pipeline(cmd, out_dir_run)

def load_latest_run_results():
    """Tự động kiểm tra các thư mục kết quả để nạp dữ liệu chạy gần nhất, hỗ trợ cả design và validation."""
    from datetime import datetime
    import pandas as pd
    
    runs_dir = os.path.join(os.getcwd(), "runs")
    if not os.path.exists(runs_dir):
        return None
        
    all_runs = []
    try:
        # Level 1 Directory: Target Species/Gene Name
        for target_name in os.listdir(runs_dir):
            target_path = os.path.join(runs_dir, target_name)
            if not os.path.isdir(target_path) or target_name.startswith("."):
                continue
            # Level 2 Directory: Run Sessions (run_YYYYMMDD_HHMMSS)
            for run_name in os.listdir(target_path):
                run_path = os.path.join(target_path, run_name)
                if not os.path.isdir(run_path) or run_name.startswith("."):
                    continue
                
                final_csv = os.path.join(run_path, "FINAL_ASSAY.csv")
                val_csv = os.path.join(run_path, "PCR_Advanced_Report.csv")
                
                # Extract timestamp
                timestamp = datetime.now()
                try:
                    parts = run_name.split("_")
                    if len(parts) >= 3:
                        ts_str = f"{parts[1]}_{parts[2]}"
                        timestamp = datetime.strptime(ts_str, "%Y%m%d_%H%M%S")
                except Exception:
                    pass
                    
                multiplex_csv = os.path.join(run_path, "MULTIPLEX_KITS.csv")
                if os.path.exists(multiplex_csv):
                    all_runs.append({
                        "path": run_path,
                        "target_name": target_name.replace("_", " ").title() + " (Multiplex)",
                        "timestamp": timestamp,
                        "type": "multiplex"
                    })
                elif os.path.exists(final_csv):
                    all_runs.append({
                        "path": run_path,
                        "target_name": target_name.replace("_", " ").title(),
                        "timestamp": timestamp,
                        "type": "design"
                    })
                elif os.path.exists(val_csv):
                    all_runs.append({
                        "path": run_path,
                        "target_name": target_name.replace("_", " ").title() + " (Validation)",
                        "timestamp": timestamp,
                        "type": "validation"
                    })
    except Exception:
        pass
                
    if not all_runs:
        # Fallback mechanism to scan old folders if runs/ is empty
        dirs_to_check = []
        if "out_dir_val" in st.session_state:
            dirs_to_check.append(st.session_state["out_dir_val"])
        if "l_out_val" in st.session_state:
            dirs_to_check.append(st.session_state["l_out_val"])
        dirs_to_check.extend(["results_auto", "results_local"])
        
        for d in dirs_to_check:
            if not d or not os.path.exists(d):
                continue
            final_csv = os.path.join(d, "FINAL_ASSAY.csv")
            multiplex_csv = os.path.join(d, "MULTIPLEX_KITS.csv")
            val_csv = os.path.join(d, "PCR_Advanced_Report.csv")
            if os.path.exists(multiplex_csv):
                all_runs.append({
                    "path": d,
                    "target_name": "Fallback Multiplex",
                    "timestamp": datetime.now(),
                    "type": "multiplex"
                })
                break
            elif os.path.exists(final_csv):
                all_runs.append({
                    "path": d,
                    "target_name": "Fallback Run",
                    "timestamp": datetime.now(),
                    "type": "design"
                })
                break
            elif os.path.exists(val_csv):
                all_runs.append({
                    "path": d,
                    "target_name": "Fallback Validation",
                    "timestamp": datetime.now(),
                    "type": "validation"
                })
                break
                
    if not all_runs:
        return None
        
    all_runs.sort(key=lambda x: x["timestamp"], reverse=True)
    latest_run = all_runs[0]
    
    results = {
        "output_dir": latest_run["path"],
        "target_name": latest_run["target_name"],
        "timestamp": latest_run["timestamp"].strftime("%Y-%m-%d %H:%M:%S"),
        "type": latest_run["type"]
    }
    
    ai_json = os.path.join(latest_run["path"], "ai_expert_report.json")
    if os.path.exists(ai_json):
        try:
            with open(ai_json, "r", encoding="utf-8") as f:
                results["expert_report"] = json.load(f)
        except Exception:
            pass
            
    if latest_run["type"] == "multiplex":
        multiplex_csv = os.path.join(latest_run["path"], "MULTIPLEX_KITS.csv")
        try:
            df = pd.read_csv(multiplex_csv)
            results["multiplex_results"] = df.to_dict(orient="records")
        except Exception:
            pass
    elif latest_run["type"] == "design":
        final_csv = os.path.join(latest_run["path"], "FINAL_ASSAY.csv")
        try:
            df = pd.read_csv(final_csv)
            columns_to_keep = ["Set_ID", "Sensitivity", "Spec", "Fwd_Primer", "Rev_Primer", "Probe_Seq", "Probe_Tm", "Target_Gene"]
            existing_cols = [c for c in columns_to_keep if c in df.columns]
            results["designed_assays"] = df[existing_cols].to_dict(orient="records")
        except Exception:
            pass
    else:
        val_csv = os.path.join(latest_run["path"], "PCR_Advanced_Report.csv")
        try:
            df = pd.read_csv(val_csv)
            results["validation_results"] = df.to_dict(orient="records")
        except Exception:
            pass
        # Load structured validation report (per-primer summary + cross-contamination)
        val_report_json = os.path.join(latest_run["path"], "validation_report.json")
        if os.path.exists(val_report_json):
            try:
                with open(val_report_json, "r", encoding="utf-8") as f:
                    val_report = json.load(f)
                results["validation_report"] = val_report
                # Inject into expert_report for AI context
                if "expert_report" not in results:
                    results["expert_report"] = {}
                results["expert_report"]["validation_summary"] = val_report

                # Prefer full cross_contamination_report.json (richer structure)
                # which has top_cross_reactive_strains, per_primer_cross_reactivity, etc.
                xc_detail_path = os.path.join(
                    latest_run["path"], "4_validation_report", "cross_contamination_report.json"
                )
                if os.path.exists(xc_detail_path):
                    with open(xc_detail_path, "r", encoding="utf-8") as f:
                        full_xc = json.load(f)
                    results["expert_report"]["cross_contamination_traceback"] = full_xc
                else:
                    # Fallback: use flattened summary from validation_report.json
                    # Remap keys to match ai_core.py expectations
                    xc_flat = val_report.get("cross_contamination", {})
                    if xc_flat:
                        results["expert_report"]["cross_contamination_traceback"] = {
                            "severity": xc_flat.get("severity", "NONE"),
                            "total_cross_reactive_strains": xc_flat.get("total_cross_reactive_strains", 0),
                            "top_cross_reactive_strains": xc_flat.get("top_offenders", []),
                            "per_primer_cross_reactivity": [],
                            "accepted_primer_offenders": [],
                            "ai_summary": xc_flat.get("ai_summary", "")
                        }
            except Exception:
                pass

            
    return results

def render_multiplex_details(run_dir, lang):
    import os
    import json
    import pandas as pd
    
    details_path = os.path.join(run_dir, "multiplex_details.json")
    ai_path = os.path.join(run_dir, "ai_expert_report.json")
    
    if not os.path.exists(details_path):
        st.warning("Không tìm thấy thông tin chi tiết tổ hợp Multiplex." if lang == "vi" else "Multiplex details file not found.")
        return
        
    try:
        with open(details_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        st.error(f"Lỗi đọc file: {e}")
        return
        
    st.markdown("### 🧬 " + ("Thermodynamic & Biophysical Analysis"))
    
    top_kits = data.get("top_kits", [])
    if not top_kits:
        st.info("Không có tổ hợp Multiplex Kit khả thi nào." if lang == "vi" else "No multiplex combinations found.")
        return
        
    kit_labels = [f"🏆 Combo #{idx+1} ({kit['verdict']} - Score: {kit['score']})" for idx, kit in enumerate(top_kits)]
    kit_tabs = st.tabs(kit_labels)
    
    for idx, kit in enumerate(top_kits):
        with kit_tabs[idx]:
            # Metric Summary
            c1, c2, c3, c4 = st.columns(4)
            score_color = "green" if kit["score"] >= 85.0 else "orange" if kit["score"] >= 70.0 else "red"
            c1.markdown(f"""
                <div style="text-align: center; background-color: #f5f5f5; border-radius: 8px; padding: 10px; border-top: 4px solid {score_color};">
                    <span style="font-size: 12px; color: gray;">{"ĐIỂM TƯƠNG THÍCH" if lang == "vi" else "COMPATIBILITY SCORE"}</span><br/>
                    <span style="font-size: 24px; font-weight: bold; color: {score_color};">{kit["score"]}/100</span>
                </div>
            """, unsafe_allow_html=True)
            
            verdict_color = "#1b5e20" if kit["verdict"] == "EXCELLENT" else "#e65100" if kit["verdict"] == "GOOD" else "#b71c1c"
            c2.markdown(f"""
                <div style="text-align: center; background-color: #f5f5f5; border-radius: 8px; padding: 10px; border-top: 4px solid {verdict_color};">
                    <span style="font-size: 12px; color: gray;">{"PHÁN QUYẾT" if lang == "vi" else "VERDICT"}</span><br/>
                    <span style="font-size: 20px; font-weight: bold; color: {verdict_color};">{kit["verdict"]}</span>
                </div>
            """, unsafe_allow_html=True)
            
            c3.markdown(f"""
                <div style="text-align: center; background-color: #f5f5f5; border-radius: 8px; padding: 10px; border-top: 4px solid #0d47a1;">
                    <span style="font-size: 12px; color: gray;">{"NHIỆT ĐỘ MỒI TRUNG BÌNH" if lang == "vi" else "MEAN PRIMER TM"}</span><br/>
                    <span style="font-size: 24px; font-weight: bold; color: #0d47a1;">{kit["mean_primer_tm"]}°C</span>
                </div>
            """, unsafe_allow_html=True)
            
            span_color = "#1b5e20" if kit["primer_tm_span"] <= 2.0 else "#e65100" if kit["primer_tm_span"] <= 4.0 else "#b71c1c"
            c4.markdown(f"""
                <div style="text-align: center; background-color: #f5f5f5; border-radius: 8px; padding: 10px; border-top: 4px solid {span_color};">
                    <span style="font-size: 12px; color: gray;">{"ĐỘ LỆCH TM MỒI" if lang == "vi" else "PRIMER TM SPAN"}</span><br/>
                    <span style="font-size: 24px; font-weight: bold; color: {span_color};">±{kit["primer_tm_span"]}°C</span>
                </div>
            """, unsafe_allow_html=True)
            
            st.markdown("<br>", unsafe_allow_html=True)
            
            # Kit Composition Table
            st.markdown("#### 🧪 " + ("Multiplex Kit Composition"))
            comp_rows = []
            for a in kit["assays"]:
                comp_rows.append({
                    "Target Name" if lang == "en" else "Tên Loài Đích": a["Target_Name"],
                    "Set ID": a["Set_ID"],
                    "Forward Primer" if lang == "en" else "Mồi Forward (5'-3')": a["Fwd_Primer"],
                    "Reverse Primer" if lang == "en" else "Mồi Reverse (5'-3')": a["Rev_Primer"],
                    "Probe Sequence" if lang == "en" else "Trình Tự Probe (5'-3')": a.get("Probe_Seq", "N/A"),
                    "Probe Tm": f"{a.get('Probe_Tm', 0.0):.1f}°C" if a.get("Probe_Seq") else "N/A",
                    "Amplicon Size (bp)" if lang == "en" else "Kích thước Sp PCR (bp)": a.get("Amplicon_Size") if a.get("Amplicon_Size") else (len(a.get("Amplicon_Sequence", "")) if a.get("Amplicon_Sequence") else "N/A")
                })
            st.dataframe(pd.DataFrame(comp_rows), use_container_width=True)
            
            # IDT Order Sheet Export
            idt_rows = []
            probe_dyes = ["/56-FAM/", "/5HEX/", "/5Cy5/", "/5TEX615/", "/5JOEN/"]
            for assay_idx, a in enumerate(kit["assays"]):
                target_base = a["Target_Name"].replace(" ", "_").replace("/", "_")
                # Fwd Primer
                idt_rows.append({
                    "Name": f"{target_base}_Fwd",
                    "Sequence": a["Fwd_Primer"],
                    "Scale": "25nm",
                    "Purification": "STD"
                })
                # Rev Primer
                idt_rows.append({
                    "Name": f"{target_base}_Rev",
                    "Sequence": a["Rev_Primer"],
                    "Scale": "25nm",
                    "Purification": "STD"
                })
                # TaqMan Probe (if exists)
                if a.get("Probe_Seq"):
                    dye = probe_dyes[assay_idx % len(probe_dyes)]
                    idt_rows.append({
                        "Name": f"{target_base}_Probe",
                        "Sequence": f"{dye}{a['Probe_Seq']}/3BHQ_1/",
                        "Scale": "100nm",
                        "Purification": "HPLC"
                    })
            
            df_idt = pd.DataFrame(idt_rows)
            import io
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                df_idt.to_excel(writer, index=False, sheet_name="IDT_Order")
            excel_data = excel_buffer.getvalue()
            
            st.download_button(
                label="📥 " + ("Tải File Excel Đặt Hàng IDT" if lang == "vi" else "Download IDT Bulk Order Excel"),
                data=excel_data,
                file_name=f"IDT_Order_Kit_{idx+1}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key=f"idt_btn_{idx}_{id(kit)}"
            )
            
            # Biophysical Warning Grid
            st.markdown("#### 🚨 " + ("Physical Biophysical Interactions & Dimer Auditing"))
            
            col_w1, col_w2 = st.columns(2)
            with col_w1:
                st.markdown("**1. Self-Dimer & Hairpins (Tương tác tự thân)**")
                if kit["self_dimer_warnings"] or kit["hairpin_warnings"]:
                    for w in kit["self_dimer_warnings"]:
                        st.warning(f"⚠️ {w}")
                    for w in kit["hairpin_warnings"]:
                        st.warning(f"🎗️ {w}")
                else:
                    st.success("✅ Không phát hiện cấu trúc Hairpin hay Self-Dimer tự thân nguy hiểm." if lang == "vi" else "✅ No critical self-folding hairpins or self-dimer structures detected.")
                    
            with col_w2:
                st.markdown("**2. Cross-Dimerization (Bắt cặp chéo giữa các mồi)**")
                if kit["cross_dimer_warnings"]:
                    for w in kit["cross_dimer_warnings"]:
                        st.error(f"🚨 {w}")
                else:
                    st.success("✅ Hoàn hảo! Không có phản ứng bám chéo chẩn đoán giữa các loài." if lang == "vi" else "✅ Perfect! No diagnostic cross-hybridization detected between targets.")
            
            # PCR Product Sizes & Gel/qPCR Separation Warnings
            st.markdown("<br>**3. " + ("Kích Thước Sản Phẩm & Phân Tách Gel" if lang == "vi" else "Product Sizing & Gel Separation Analysis") + "**", unsafe_allow_html=True)
            if kit.get("size_warnings"):
                for w in kit["size_warnings"]:
                    if "🚨" in w:
                        st.error(w)
                    else:
                        st.warning(w)
            else:
                st.success("✅ Kích thước các sản phẩm PCR hoàn hảo và tương thích tuyệt đối với phương pháp đã chọn." if lang == "vi" else "✅ PCR product sizes are perfect and fully compatible with the selected assay format.")
            
    # AI Clinical Report Panel
    if os.path.exists(ai_path):
        try:
            with open(ai_path, "r", encoding="utf-8") as f:
                rep = json.load(f)
            
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown(f"""
                <div style="background-color: #e3f2fd; border-left: 5px solid #0d47a1; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05);">
                    <span style="font-weight: bold; color: #0d47a1; font-size: 18px;">🤖 AI Expert Clinical Advice & qPCR Wet-lab Recommendations</span>
                </div>
            """, unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)
            
            col_ai1, col_ai2 = st.columns([1, 2])
            with col_ai1:
                st.metric(label="🏆 " + ("BỘ KIT ĐỀ XUẤT" if lang == "vi" else "RECOMMENDED KIT"), value=rep.get("best_assay_name", "N/A"))
                st.metric(label="📊 " + ("PHÊ DUYỆT LÂM SÀNG" if lang == "vi" else "CLINICAL VERDICT"), value=rep.get("overall_verdict", "N/A"))
                if 'confidence_score' in rep:
                    st.metric(label="💡 " + ("ĐỘ TỰ TIN (AI)" if lang == "vi" else "AI CONFIDENCE"), value=f"{rep['confidence_score']}%")
            with col_ai2:
                st.info(f"**🌡️ Cân Bằng Tm & GC:** {rep.get('tm_analysis', 'N/A')}")
                st.success(f"**🎯 Phân Tích Độ Đặc Hiệu:** {rep.get('specificity_sensitivity_balance', 'N/A')}")
                if 'decision_justification' in rep:
                    st.warning(f"**Lập luận thuật toán (Algorithm Justification):** {rep['decision_justification']}")
                
            st.markdown(f"**💡 Hướng dẫn Tối ưu Hóa Wet-Lab:**\n{rep.get('clinical_recommendation', 'N/A')}")
            
            if rep.get("structural_risks"):
                for w in rep.get("structural_risks"):
                    st.warning(f"⚠️ {w}")
        except Exception:
            pass

# --- 7. TABS CONTENT ---

tab_history, tab_local, tab_auto, tab_chat = st.tabs([
    tr("tab_history_title"),
    tr("tab_local_title"),
    tr("tab_auto_title"),
    tr("tab_chat_title")
])

with tab_chat:
    st.info(tr("chat_suggest"))
    import re
    
    # Automatically load the latest run result if available
    run_results = load_latest_run_results()
    
    if run_results:
        is_val = run_results.get("type") == "validation"
        is_multi = run_results.get("type") == "multiplex"
        if is_multi:
            title_text = "🧬 Trạng thái Chuyên gia AI (AI Expert): Kết Quả Multiplex Đã Nạp!" if lang == "vi" else "🧬 AI Expert Status: Multiplex Results Loaded!"
            desc_text = f"Phát hiện kết quả thiết kế Multiplex thành công tại thư mục <code>{run_results['output_dir']}</code>. AI đã tự động nạp dữ liệu phân tích tương thích tổ hợp mồi/probe và sẵn sàng tư vấn chuyên sâu!" if lang == "vi" else f"Successfully loaded multiplex assay design from <code>{run_results['output_dir']}</code>. AI is ready to provide clinical multiplex wet-lab advice!"
        elif is_val:
            title_text = tr("validation_loaded")
            desc_text = tr("validation_desc").format(dir=run_results['output_dir'])
        else:
            title_text = tr("results_loaded")
            desc_text = tr("results_desc").format(dir=run_results['output_dir'])
        
        st.markdown(
            f"""
            <div style="background-color: #f0fdf4; border: 1px solid #bbf7d0; border-left: 5px solid #16a34a; padding: 16px; border-radius: 12px; margin-bottom: 20px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.03);">
                <span style="font-family: 'Outfit', sans-serif; font-weight: 600; color: #166534; font-size: 16px;">{title_text}</span><br/>
                <span style="font-family: 'Outfit', sans-serif; font-size: 14px; color: #334155; margin-top: 5px; display: inline-block; line-height: 1.4;">{desc_text}</span>
            </div>
            """, 
            unsafe_allow_html=True
        )
        
        # Display corresponding results table
        if is_multi and "multiplex_results" in run_results:
            pass # Skip displaying multiplex table on chat/running page, only show in dashboard
        elif not is_val and "designed_assays" in run_results:
            with st.expander(tr("quick_assays_title"), expanded=True):
                import pandas as pd
                df_assays = pd.DataFrame(run_results["designed_assays"])
                st.dataframe(df_assays, use_container_width=True)
        elif is_val and "validation_results" in run_results:
            with st.expander(tr("quick_validation_title"), expanded=True):
                import pandas as pd
                df_val = pd.DataFrame(run_results["validation_results"])
                cols = [c for c in ["id", "group", "marker", "status", "mismatch", "amplicon_size", "pair_type", "copy_number", "identity"] if c in df_val.columns]
                st.dataframe(df_val[cols], use_container_width=True)
        
        col_reset, _ = st.columns([1.5, 3])
        with col_reset:
            if st.button(tr("btn_reset"), use_container_width=True, key="btn_reset_assays"):
                try:
                    final_csv = os.path.join(run_results['output_dir'], "FINAL_ASSAY.csv")
                    val_csv = os.path.join(run_results['output_dir'], "PCR_Advanced_Report.csv")
                    ai_json = os.path.join(run_results['output_dir'], "ai_expert_report.json")
                    if os.path.exists(final_csv): os.remove(final_csv)
                    if os.path.exists(val_csv): os.remove(val_csv)
                    if os.path.exists(ai_json): os.remove(ai_json)
                    st.rerun()
                except Exception as e:
                    st.error(f"{tr('error_reset')}{e}")
    
    chat_container = st.container(height=500)
    with chat_container:
        # Display chat history
        for msg in st.session_state["chat_history"]:
            with st.chat_message(msg["role"]):
                st.markdown(msg["content"])
            
    # Display friendly recommended configuration
    json_placeholder = st.empty()
    run_clicked = False
    if st.session_state["copilot_json"]:
        with json_placeholder.container():
            cfg = st.session_state["copilot_json"]
            st.info(tr("ai_setup_config"))
            
            if cfg.get("action") == "propose_multiplex":
                c_t, c_b = st.columns(2)
                with c_t:
                    st.markdown("🎯 **" + ("Các loài Target (Multiplex Targets)" if lang == "vi" else "Multiplex Targets") + "**")
                    for t in cfg.get("targets", []):
                        st.write(f"- 🧬 *{t['query']}* ({t.get('count', 500)} {'strains' if lang == 'en' else 'mẫu'})")
                with c_b:
                    st.markdown("🛡️ **" + ("Bộ nền loại trừ chung (Shared Background)" if lang == "vi" else "Shared Background") + "**")
                    for b in cfg.get("background", []):
                        st.write(f"- 🦠 *{b['query']}* ({b.get('count', 50)} {'strains' if lang == 'en' else 'mẫu'})")
            else:
                c_t, c_b = st.columns(2)
                with c_t:
                    st.markdown(tr("target_inclusion"))
                    for t in cfg.get("target", []):
                        size_str = f" >= {t['size']} Mb" if t.get('size', 0) > 0 else ""
                        st.write(f"- 🧬 *{t['query']}* ({t.get('count', 50)} {'strains' if lang == 'en' else 'mẫu'}{size_str})")
                with c_b:
                    st.markdown(tr("bg_exclusion"))
                    for b in cfg.get("background", []):
                        size_str = f" >= {b['size']} Mb" if b.get('size', 0) > 0 else ""
                        st.write(f"- 🦠 *{b['query']}* ({b.get('count', 10)} {'strains' if lang == 'en' else 'mẫu'}{size_str})")
                        
                if cfg.get("action") == "propose_validation":
                    st.markdown(tr("primers_to_validate_label"))
                    for p in cfg.get("primers", []):
                        st.write(f"- 🏷️ **{p['name']}**: `5'-{p['fwd']}-3'` (Fwd) / `5'-{p['rev']}-3'` (Rev)")
                    
            st.markdown("---")
            confirm_tip = (
                "👉 *Do you want to adjust or add any species? If not, click the button below or type 'run pipeline' to start.*"
                if lang == "en" else
                "👉 *Bạn có muốn thay đổi hay thêm loài nào không? Nếu không, chỉ cần nhấn nút bên dưới hoặc gõ 'chạy pipeline' để bắt đầu.*"
            )
            st.markdown(confirm_tip)
            
            if st.button(tr("confirm_run"), type="primary", key="btn_run_copilot"):
                run_clicked = True
            else:
                run_clicked = False
                
    if run_clicked:
        if not st.session_state["email_val"] or "@" not in st.session_state["email_val"]:
            st.error(tr("error_email"))
        else:
            json_placeholder.empty() # Clear UI immediately so user knows it started
            cfg = st.session_state["copilot_json"]
            execute_proposal_pipeline(cfg)
            st.session_state["copilot_json"] = None # Reset
            
            with chat_container:
                if st.session_state["chat_history"]:
                    msg = st.session_state["chat_history"][-1]
                    if msg["role"] == "assistant" and "Tôi đã đánh giá xong" in msg["content"]:
                        with st.chat_message(msg["role"]):
                            st.markdown(msg["content"])

    # Chat input box
    if prompt := st.chat_input(tr("chat_input_placeholder")):
        st.session_state["chat_history"].append({"role": "user", "content": prompt})
        
        with chat_container:
            with st.chat_message("user"):
                st.markdown(prompt)
                
            # Automatically trigger pipeline if user types a command
            is_trigger_command = prompt.strip().lower() in (
                ["chạy pipeline", "chay pipeline", "run pipeline", "chạy đi", "ok chạy", "ok run", "start pipeline"]
            )
            if is_trigger_command and st.session_state["copilot_json"]:
                if not st.session_state["email_val"] or "@" not in st.session_state["email_val"]:
                    st.error(tr("error_email"))
                else:
                    json_placeholder.empty() # Clear UI
                    cfg = st.session_state["copilot_json"]
                    execute_proposal_pipeline(cfg)
                    st.session_state["copilot_json"] = None
                    
                    if st.session_state["chat_history"]:
                        msg = st.session_state["chat_history"][-1]
                        if msg["role"] == "assistant" and "đánh giá xong" in msg["content"].lower():
                            with st.chat_message(msg["role"]):
                                st.markdown(msg["content"])
            else:
                with st.chat_message("assistant"):
                    with st.spinner("Thinking..." if lang == "en" else "Đang suy nghĩ..."):
                        try:
                            from rational_design.ai_core import AIExpertAgent, LocalBackend
                            from rational_design.analytics import ResultAnalysisEngine
                            
                            backend = LocalBackend(base_url=ai_base_url, model_name=ai_model)
                            agent = AIExpertAgent(backend)
                            
                            # Index all historical runs to load context for AI
                            engine = ResultAnalysisEngine("runs")
                            past_runs = engine.scan_historical_runs()
                            history_summary = engine.generate_ai_summary(past_runs)
                            
                            copilot_context = {
                                "latest_completed_run": run_results,
                                "historical_database_summary": history_summary
                            }
                            
                            history_for_api = [{"role": m["role"], "content": m["content"]} for m in st.session_state["chat_history"]]
                            
                            # Stream response and pass copilot_context
                            response_stream = agent.chat_stream(history_for_api, expert_report=copilot_context, language=lang)
                            response_text = ""
                            
                            if hasattr(st, "write_stream"):
                                response_text = st.write_stream(response_stream)
                            else:
                                placeholder = st.empty()
                                for chunk in response_stream:
                                    response_text += chunk
                                    placeholder.markdown(response_text + "▌")
                                placeholder.markdown(response_text)
                            
                            json_match = re.search(r'```json\s*(\{.*?\})\s*```', response_text, re.DOTALL)
                            if not json_match:
                                json_match = re.search(r'(\{[\s\S]*"action":\s*"(propose_design|propose_local_design|propose_validation|propose_multiplex)"[\s\S]*\})', response_text)
        
                            display_text = response_text
                            parsed_json = None
                            if json_match:
                                try:
                                    parsed_json = json.loads(json_match.group(1))
                                    display_text = response_text.replace(json_match.group(0), "")
                                except Exception:
                                    pass
                                    
                            st.session_state["chat_history"].append({"role": "assistant", "content": display_text})
                            
                            if parsed_json and parsed_json.get("action") in ["propose_design", "propose_local_design", "propose_validation", "propose_multiplex"]:
                                if parsed_json.get("run_immediately") is True:
                                    if not st.session_state["email_val"] or "@" not in st.session_state["email_val"]:
                                        st.error(tr("error_email"))
                                    else:
                                        json_placeholder.empty()
                                        execute_proposal_pipeline(parsed_json)
                                else:
                                    st.session_state["copilot_json"] = parsed_json
                                    st.session_state["chat_history"].append({
                                        "role": "assistant", 
                                        "content": tr("chat_history_finish")
                                    })
                                    st.rerun()
                        except Exception as e:
                            st.error(f"{'AI connection error' if lang == 'en' else 'Lỗi kết nối AI'}: {e}")
with tab_auto:
    st.subheader(tr("auto_header"))
    col1, col2 = st.columns(2)
    email = col1.text_input(tr("ncbi_email_disabled"), placeholder=tr("sidebar_placeholder"), disabled=True)
    project_name = col2.text_input(tr("project_name"), key="auto_proj_val")
    
    is_multiplex = st.checkbox("Thiết Kế Multiplex PCR (Đa Tác Nhân)" if lang == "vi" else "Multiplex PCR Mode (Multi-target)", value=st.session_state["is_multiplex"], key="is_multiplex_cb")
    st.session_state["is_multiplex"] = is_multiplex
    
    if is_multiplex:
        st.session_state["multiplex_type"] = st.radio(
            "Phương pháp Multiplex PCR (Multiplex Assay Type)" if lang == "vi" else "Multiplex Assay Type",
            options=["qPCR", "Conventional"],
            format_func=lambda x: ("TaqMan qPCR (Probe-based)" if x == "qPCR" else "Conventional PCR (Gel-based)") if lang == "vi" else ("TaqMan qPCR (Probe-based)" if x == "qPCR" else "Conventional PCR (Gel-based)"),
            index=0 if st.session_state.get("multiplex_type", "qPCR") == "qPCR" else 1,
            horizontal=True
        )
    
    st.markdown("---")
    if is_multiplex:
        st.markdown("#### 🎯 Individual Multiplex Targets")
        for i, item in enumerate(st.session_state["multiplex_target_list"]):
            c1, c2, c3, c4 = st.columns([5, 2, 2, 1])
            with c1:
                item["query"] = st.text_input(tr("target_query_num").format(num=i+1), value=item["query"], key=f"t_multi_q_{i}_{st.session_state.reset_id}")
            with c2:
                item["size"] = st.number_input(tr("min_mb_num").format(num=i+1), value=item["size"], step=0.1, key=f"t_multi_s_{i}_{st.session_state.reset_id}")
            with c3:
                item["count"] = st.number_input(tr("max_count_num").format(num=i+1), value=item.get("count", 500), min_value=0, key=f"t_multi_c_{i}_{st.session_state.reset_id}")
            with c4:
                st.markdown("<div style='height: 28px;'></div>", unsafe_allow_html=True)
                if st.button("🗑️", key=f"t_multi_del_{i}_{st.session_state.reset_id}"):
                    if len(st.session_state["multiplex_target_list"]) > 2:
                        st.session_state["multiplex_target_list"].pop(i)
                        st.rerun()
                        
        if st.button(tr("add_multiplex_target"), key="add_multiplex_target"):
            st.session_state["multiplex_target_list"].append({"query": "", "size": 0.0, "count": 500})
            st.rerun()
    else:
        st.markdown(tr("target_group_header"))
        for i, item in enumerate(st.session_state["target_list"]):
            c1, c2, c3, c4 = st.columns([5, 2, 2, 1])
            with c1:
                item["query"] = st.text_input(tr("target_query_num").format(num=i+1), value=item["query"], key=f"t_q_{i}_{st.session_state.reset_id}")
            with c2:
                item["size"] = st.number_input(tr("min_mb_num").format(num=i+1), value=item["size"], step=0.1, key=f"t_s_{i}_{st.session_state.reset_id}")
            with c3:
                item["count"] = st.number_input(tr("max_count_num").format(num=i+1), value=item.get("count", 500), min_value=0, key=f"t_c_{i}_{st.session_state.reset_id}")
            with c4:
                st.markdown("<div style='height: 28px;'></div>", unsafe_allow_html=True)
                if st.button("🗑️", key=f"t_del_{i}_{st.session_state.reset_id}"):
                    if len(st.session_state["target_list"]) > 1:
                        st.session_state["target_list"].pop(i)
                        st.rerun()
     
        if st.button(tr("add_target"), key="add_target"):
            st.session_state["target_list"].append({"query": "", "size": 0.0, "count": 500})
            st.rerun()
 
    st.markdown("---")
    st.markdown(tr("bg_group_header"))
    for i, item in enumerate(st.session_state["bg_list"]):
        c1, c2, c3, c4 = st.columns([5, 2, 2, 1])
        with c1:
            item["query"] = st.text_input(tr("bg_query_num").format(num=i+1), value=item["query"], key=f"b_q_{i}_{st.session_state.reset_id}")
        with c2:
            item["size"] = st.number_input(tr("min_mb_num").format(num=i+1), value=item["size"], step=0.1, key=f"b_s_{i}_{st.session_state.reset_id}")
        with c3:
            item["count"] = st.number_input(tr("max_count_num").format(num=i+1), value=item.get("count", 50), min_value=0, key=f"b_c_{i}_{st.session_state.reset_id}")
        with c4:
            st.markdown("<div style='height: 28px;'></div>", unsafe_allow_html=True)
            if st.button("🗑️", key=f"b_del_{i}_{st.session_state.reset_id}"):
                if len(st.session_state["bg_list"]) > 1:
                    st.session_state["bg_list"].pop(i)
                    st.rerun()
                    
    if st.button(tr("add_bg"), key="add_bg"):
        st.session_state["bg_list"].append({"query": "", "size": 0.0, "count": 50})
        st.rerun()
    
    st.divider()
    if st.button(tr("btn_start_auto"), type="primary"):
        if not st.session_state["email_val"] or "@" not in st.session_state["email_val"]:
            st.error(tr("error_email_valid"))
        else:
            param_file = save_params()
            if is_multiplex:
                import re
                from datetime import datetime
                
                valid_targets = [t for t in st.session_state["multiplex_target_list"] if t["query"]]
                if len(valid_targets) < 2:
                    st.error("Cần nhập tối thiểu 2 loài Target để thiết kế Multiplex." if lang == "vi" else "At least 2 target queries are required for multiplex design.")
                else:
                    target_names = [t["query"] for t in valid_targets]
                    proj_name = project_name.strip() if project_name.strip() else "multiplex_" + "_".join([re.sub(r'[^a-zA-Z0-9]', '', n.replace(" ", "_").lower()) for n in target_names])
                    if len(proj_name) > 60:
                        proj_name = proj_name[:60]
                        
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    out_dir_run = os.path.join(os.getcwd(), "runs", proj_name, f"run_{timestamp}")
                    
                    cfg = {
                        "action": "propose_multiplex",
                        "targets": [{"query": t["query"], "count": t.get("count", 500), "size": t.get("size", 0.0), "type": "genome"} for t in valid_targets],
                        "background": [{"query": b["query"], "count": b.get("count", 50), "size": b.get("size", 0.0), "type": "genome"} for b in st.session_state["bg_list"] if b["query"]]
                    }
                    execute_proposal_pipeline(cfg, out_dir_run)
            else:
                t_conf = {f"t{idx+1}": [it["query"], it["size"], it["count"]] for idx, it in enumerate(st.session_state["target_list"]) if it["query"]}
                b_conf = {f"b{idx+1}": [it["query"], it["size"], it["count"]] for idx, it in enumerate(st.session_state["bg_list"]) if it["query"]}
                
                with open("config/t_conf.json", "w") as f: json.dump(t_conf, f)
                with open("config/b_conf.json", "w") as f: json.dump(b_conf, f)
                
                target_species_name = st.session_state["target_list"][0]["query"] if st.session_state["target_list"] else "unknown"
                out_dir_run = generate_run_output_dir(target_species_name)
                
                cmd = f'"{sys.executable}" -m rational_design.cli pipeline --out "{out_dir_run}" --email "{st.session_state["email_val"]}" --target_config "config/t_conf.json" --bg_config "config/b_conf.json" --params "{param_file}" --language "{lang}"'
                if enable_ai:
                    cmd += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
                run_pipeline(cmd, out_dir_run)

with tab_local:
    st.subheader(tr("local_header"))
    l_proj = st.text_input(tr("project_name"), key="local_proj_val", value="T_Marneffei_Run")
    
    is_local_multiplex = st.checkbox(
        "Thiết Kế Multiplex PCR (Đa Tác Nhân Local)" if lang == "vi" else "Local Multiplex PCR Mode (Multi-target)", 
        value=st.session_state.get("is_local_multiplex", False), 
        key="is_local_multiplex_cb"
    )
    st.session_state["is_local_multiplex"] = is_local_multiplex
    
    if is_local_multiplex:
        st.session_state["multiplex_type"] = st.radio(
            "Phương pháp Multiplex PCR (Multiplex Assay Type)" if lang == "vi" else "Multiplex Assay Type",
            options=["qPCR", "Conventional"],
            format_func=lambda x: ("TaqMan qPCR (Probe-based)" if x == "qPCR" else "Conventional PCR (Gel-based)") if lang == "vi" else ("TaqMan qPCR (Probe-based)" if x == "qPCR" else "Conventional PCR (Gel-based)"),
            index=0 if st.session_state.get("multiplex_type", "qPCR") == "qPCR" else 1,
            horizontal=True,
            key="local_multiplex_type_radio"
        )
        st.markdown("#### 🎯 Local Target Folders (Inclusion)")
        if "local_target_list" not in st.session_state:
            st.session_state["local_target_list"] = ["", ""]
            
        for i in range(len(st.session_state["local_target_list"])):
            c1, c2, c3 = st.columns([5, 1, 1])
            with c1:
                st.session_state["local_target_list"][i] = st.text_input(
                    f"Đường dẫn Target Local #{i+1}" if lang == "vi" else f"Local Target Folder #{i+1}",
                    value=st.session_state["local_target_list"][i],
                    key=f"l_t_{i}_{st.session_state.reset_id}"
                )
            with c2:
                st.markdown("<div style='height: 28px;'></div>", unsafe_allow_html=True)
                if HAS_TK:
                    st.button(tr("browse"), key=f"btn_l_t_{i}", on_click=select_folder_callback, args=(f"l_t_{i}_temp_key",))
                    if st.session_state.get(f"l_t_{i}_temp_key"):
                        st.session_state["local_target_list"][i] = st.session_state.get(f"l_t_{i}_temp_key")
                        st.session_state[f"l_t_{i}_temp_key"] = None # Reset
                        st.rerun()
            with c3:
                st.markdown("<div style='height: 28px;'></div>", unsafe_allow_html=True)
                if st.button("🗑️", key=f"l_t_del_{i}_{st.session_state.reset_id}"):
                    if len(st.session_state["local_target_list"]) > 2:
                        st.session_state["local_target_list"].pop(i)
                        st.rerun()
                        
        if st.button("➕ Thêm đường dẫn Target Local" if lang == "vi" else "➕ Add Local Target Folder", key="add_local_target_btn"):
            st.session_state["local_target_list"].append("")
            st.rerun()
    else:
        col_t1, col_t2 = st.columns([4, 1])
        path_t = col_t1.text_input(tr("local_target_path"), key="path_t_val")
        if HAS_TK: col_t2.button(tr("browse"), key="btn_t", on_click=select_folder_callback, args=("path_t_val",))

    col_b1, col_b2 = st.columns([4, 1])
    path_b = col_b1.text_input(tr("local_bg_path"), key="path_b_val")
    if HAS_TK: col_b2.button(tr("browse"), key="btn_b", on_click=select_folder_callback, args=("path_b_val",))
    
    col_o1, col_o2 = st.columns([4, 1])
    l_out = col_o1.text_input(tr("local_out_path"), key="l_out_val", placeholder=tr("runs_placeholder"))
    if HAS_TK: col_o2.button(tr("browse"), key="btn_o", on_click=select_folder_callback, args=("l_out_val",))
    
    st.divider()
    if st.button(tr("btn_start_local"), type="primary"):
        if is_local_multiplex:
            import re
            import shutil
            from datetime import datetime
            
            valid_targets = [p.strip() for p in st.session_state["local_target_list"] if p.strip()]
            invalid_paths = [p for p in valid_targets if not os.path.exists(p)]
            
            if len(valid_targets) < 2:
                st.error("Cần nhập tối thiểu 2 thư mục Target để thiết kế Multiplex Local." if lang == "vi" else "At least 2 local target directories are required for local multiplex.")
            elif invalid_paths:
                st.error(f"{tr('invalid_target')}: {', '.join(invalid_paths)}")
            else:
                target_names = [os.path.basename(p) for p in valid_targets]
                proj_name = l_proj.strip() if l_proj.strip() else "local_multiplex_" + "_".join([re.sub(r'[^a-zA-Z0-9]', '', n.replace(" ", "_").lower()) for n in target_names])
                if len(proj_name) > 60:
                    proj_name = proj_name[:60]
                    
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                out_dir_run = os.path.join(os.getcwd(), "runs", proj_name, f"run_{timestamp}")
                os.makedirs(out_dir_run, exist_ok=True)
                
                param_file = save_params()
                shared_bg_dir = path_b.strip()
                
                cmd_parts = []
                target_folders = []
                shared_context_file = os.path.join(out_dir_run, "shared_primer_context.json")
                
                for idx, t_path in enumerate(valid_targets):
                    target_folder = os.path.join(out_dir_run, f"t{idx+1}")
                    target_folders.append(target_folder)
                    os.makedirs(target_folder, exist_ok=True)
                    
                    # Proactively save target metadata for multiplex analyzer
                    try:
                        t_name = os.path.basename(t_path)
                        with open(os.path.join(target_folder, "metadata.json"), "w", encoding="utf-8") as mf:
                            json.dump({"target_name": t_name}, mf, ensure_ascii=False, indent=4)
                    except Exception:
                        pass
                    
                    virtual_bg_dir = os.path.join(out_dir_run, "0_raw_data", f"bg_for_t{idx+1}")
                    os.makedirs(virtual_bg_dir, exist_ok=True)
                    
                    # Copy shared backgrounds
                    if shared_bg_dir and os.path.exists(shared_bg_dir):
                        for item in os.listdir(shared_bg_dir):
                            s_path = os.path.join(shared_bg_dir, item)
                            if os.path.isfile(s_path):
                                shutil.copy2(s_path, os.path.join(virtual_bg_dir, item))
                                
                    # Copy other targets
                    for jdx, other_t_path in enumerate(valid_targets):
                        if jdx != idx:
                            if os.path.exists(other_t_path):
                                for item in os.listdir(other_t_path):
                                    s_path = os.path.join(other_t_path, item)
                                    if os.path.isfile(s_path):
                                        shutil.copy2(s_path, os.path.join(virtual_bg_dir, f"target_bg_t{jdx+1}_{item}"))
                                        
                    cmd_part = f'"{sys.executable}" -m rational_design.cli pipeline --out "{target_folder}" --local_target "{t_path}" --local_bg "{virtual_bg_dir}" --params "{param_file}" --language "{lang}"'
                    if enable_ai:
                        cmd_part += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
                    cmd_part += f' --shared_context "{shared_context_file}"'
                    cmd_parts.append(cmd_part)
                    
                assay_type = st.session_state.get("multiplex_type", "qPCR")
                folders_str = " ".join([f'"{f}"' for f in target_folders])
                cmd_multi = f'"{sys.executable}" -m rational_design.cli multiplex_analyze --folders {folders_str} --out "{out_dir_run}" --language "{lang}" --assay_type "{assay_type}"'
                if enable_ai:
                    cmd_multi += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
                cmd_parts.append(cmd_multi)
                
                cmd = " && ".join(cmd_parts)
                run_pipeline(cmd, out_dir_run)
        else:
            path_t = st.session_state.get("path_t_val", "")
            if not path_t or not os.path.exists(path_t):
                st.error(tr("invalid_target"))
            else:
                param_file = save_params()
                target_species_name = os.path.basename(path_t)
                out_dir_run = generate_run_output_dir(target_species_name)
                
                cmd = f'"{sys.executable}" -m rational_design.cli pipeline --out "{out_dir_run}" --local_target "{path_t}" --local_bg "{path_b}" --params "{param_file}" --language "{lang}"'
                if enable_ai:
                    cmd += f' --ai_base_url "{ai_base_url}" --ai_model "{ai_model}"'
                run_pipeline(cmd, out_dir_run)

with tab_history:
    st.subheader(tr("history_header"))
    st.markdown(tr("history_subheader"))
    
    from rational_design.analytics import ResultAnalysisEngine
    engine = ResultAnalysisEngine("runs")
    past_runs = engine.scan_historical_runs()
    
    if not past_runs:
        st.info(tr("no_history"))
    else:
        # Display overview statistics
        c1, c2, c3 = st.columns(3)
        c1.metric(tr("total_targets"), len(set(r["target_name"] for r in past_runs)))
        c2.metric(tr("total_runs"), len(past_runs))
        c3.metric(tr("latest_run"), past_runs[0]["target_name"])
        
        st.divider()
        
        # History table as Pandas DataFrame
        import pandas as pd
        run_data = []
        for r in past_runs:
            run_data.append({
                tr("df_species"): r["target_name"],
                tr("df_time"): r["timestamp"].strftime("%Y-%m-%d %H:%M:%S"),
                tr("df_best"): r["best_assay"],
                tr("df_sensitivity"): r["sensitivity"],
                tr("df_specificity"): r["specificity"],
                tr("df_target_gene"): r["target_gene"],
                tr("df_candidates"): r["total_candidates"],
                tr("df_folder"): r["run_folder_name"]
            })
        df_runs = pd.DataFrame(run_data)
        st.dataframe(df_runs, use_container_width=True)
        
        st.divider()
        st.markdown(tr("detail_explorer"))
        
        # Run selection
        run_options = {f"{r['target_name']} ({r['timestamp'].strftime('%Y-%m-%d %H:%M:%S')})": r for r in past_runs}
        selected_run_label = st.selectbox(tr("select_run"), list(run_options.keys()))
        
        if selected_run_label:
            selected_run = run_options[selected_run_label]
            details = engine.get_run_details(selected_run["path"])
            
            st.markdown(tr("detail_run_for").format(target=selected_run['target_name']))
            st.markdown(tr("detail_path").format(path=selected_run['path']))
            
            if details.get("is_multiplex"):
                multiplex_csv = os.path.join(selected_run["path"], "MULTIPLEX_KITS.csv")
                if os.path.exists(multiplex_csv):
                    with st.expander(tr("quick_multiplex_assays_title"), expanded=True):
                        df_multi = pd.read_csv(multiplex_csv)
                        st.dataframe(df_multi, use_container_width=True)
                render_multiplex_details(selected_run["path"], lang)
            else:
                # AI Expert evaluation
                if details.get("report"):
                    rep = details["report"]
                    st.info(f"{tr('ai_recommendation')}{rep.get('clinical_recommendation', 'N/A')}")
                    col_sub1, col_sub2 = st.columns(2)
                    col_sub1.markdown(f"{tr('best_primer_pair')}{rep.get('best_assay_name', 'N/A')}")
                    col_sub2.markdown(f"{tr('overall_evaluation')}{rep.get('overall_verdict', 'N/A')}")
                    
                # Display list of Assays
                if details.get("assays"):
                    df_assays = pd.DataFrame(details["assays"])
                    cols_to_show = [c for c in ["Set_ID", "Sensitivity", "Spec", "Fwd_Primer", "Rev_Primer", "Probe_Seq", "Probe_Tm", "Target_Gene"] if c in df_assays.columns]
                    st.dataframe(df_assays[cols_to_show], use_container_width=True)
                else:
                    st.warning(tr("no_assays_found"))
