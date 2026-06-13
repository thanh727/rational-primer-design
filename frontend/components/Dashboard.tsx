"use client";

import Link from "next/link";
import { FormEvent, ReactNode, useEffect, useMemo, useRef, useState } from "react";

type Lang = "vi" | "en";
type WorkspaceView = "dashboard" | "local" | "auto" | "ai" | "validate" | "multiplex" | "about";
type DataMode = "local" | "auto";
type MultiplexMode = "local" | "auto" | "analyze";
type BrowseKind = "file" | "directory" | "any";

type Job = {
  id: string;
  status: string;
  source?: string;
  output_dir: string;
  created_at: number;
  updated_at: number;
  command: string[];
  return_code?: number | null;
};

type QueryItem = {
  query: string;
  size: number;
  count: number;
  type: "genome" | "gene";
};

type PrimerPair = {
  name: string;
  fwd: string;
  rev: string;
  probe?: string;
};

type ChatMessage = {
  role: "user" | "assistant";
  content: string;
};

type CsvPreview = {
  exists: boolean;
  columns: string[];
  rows: Record<string, unknown>[];
  row_count?: number;
  error?: string;
};

type ResultPayload = {
  job: Job;
  summary: Record<string, unknown> | null;
  audit: Record<string, unknown> | null;
  cross_contamination: Record<string, unknown> | null;
  ai_report: Record<string, unknown> | null;
  final_assays: CsvPreview;
  candidate_summary: CsvPreview;
  known_primer_validation: CsvPreview;
  known_primer_summary: unknown;
  multiplex_kits: CsvPreview;
  multiplex_details: Record<string, unknown> | null;
  files: Array<{ name: string; path: string; size_bytes: number }>;
};

type HistoryPayload = {
  web_jobs: Job[];
};

type ChatResponse = {
  reply: string;
  raw_reply: string;
  proposal: Record<string, unknown> | null;
  job: Job | null;
  blocked_reason: string | null;
};

type Parameters = Record<string, string | number | boolean>;

type OfflineModelsResponse = {
  base_url: string;
  provider: string;
  available: boolean;
  models: string[];
  error: string | null;
};

type FileBrowserEntry = {
  name: string;
  path: string;
  is_dir: boolean;
  size_bytes: number | null;
};

type FileBrowserResponse = {
  path: string;
  parent: string | null;
  entries: FileBrowserEntry[];
  roots: string[];
};

type FileBrowserState = {
  open: boolean;
  path: string;
  kind: BrowseKind;
  onSelect: (path: string) => void;
};

type StatusPayload = {
  status: string;
  run_root: string;
  file_browser_roots: string[];
  job_counts: Record<string, number>;
  job_count: number;
  run_root_size_bytes: number;
  disk: { total: number; used: number; free: number };
  ai: OfflineModelsResponse;
};

function getApiBase(): string {
  if (typeof window !== "undefined" && (window as any).electronAPI) {
    return "http://127.0.0.1:8000";
  }
  return process.env.NEXT_PUBLIC_API_BASE_URL || "http://127.0.0.1:8000";
}
const API_BASE = getApiBase();

const I18N: Record<Lang, Record<string, string>> = {
  vi: {
    app: "Rational Primer Design",
    subtitle: "Không gian thiết kế assay qPCR",
    sharedConfig: "Cấu hình dùng chung",
    aiConfig: "Chuyên gia AI",
    dashboard: "Bảng điều khiển",
    dashboardTitle: "Bảng điều khiển",
    local: "Tệp cục bộ",
    auto: "Từ khóa tự động",
    ai: "Trò chuyện AI",
    validate: "Kiểm chứng mồi",
    multiplex: "Đa mồi",
    about: "Giới thiệu",
    aboutTitle: "Giới thiệu Rational Primer Design",
    aboutIntro:
      "Rational Primer Design là không gian thiết kế, kiểm chứng và phân tích assay PCR/qPCR tích hợp, kết hợp pipeline tin sinh học offline với AI Expert chuyên sâu về sinh học phân tử.",
    aboutDesignTitle: "Thiết kế assay có kiểm soát",
    aboutDesignText:
      "Hệ thống hỗ trợ thiết kế từ FASTA cục bộ, tự động tải dữ liệu bằng từ khóa NCBI, kiểm chứng bộ mồi đã biết và xây dựng panel multiplex. Các tham số lõi được gom trong sidebar để dùng nhất quán trên toàn bộ tác vụ.",
    aboutAiTitle: "AI Expert cho PCR/qPCR",
    aboutAiText:
      "AI không chỉ tạo cấu hình chạy pipeline. Trợ lý có thể tư vấn chọn target/background, đánh giá primer/probe, diễn giải kết quả in-silico, nhận diện rủi ro dimer/hairpin/cross-reactivity và đề xuất kế hoạch tối ưu wet-lab.",
    aboutLocalTitle: "Ưu tiên dữ liệu và model nội bộ",
    aboutLocalText:
      "Package có thể chạy với dữ liệu local, duyệt thư mục trực tiếp từ giao diện, tự nhận diện model AI offline tương thích OpenAI/Ollama và lưu lịch sử job để theo dõi kết quả thiết kế.",
    aboutWorkflowTitle: "Luồng làm việc chính",
    aboutWorkflowDesign: "Thiết kế singleplex hoặc multiplex từ genome target/background.",
    aboutWorkflowValidate: "Kiểm chứng primer đã biết bằng mô phỏng PCR nâng cao.",
    aboutWorkflowReview: "Theo dõi log, bảng kết quả, báo cáo AI và file đầu ra trên dashboard.",
    history: "Lịch sử",
    monitor: "Theo dõi",
    results: "Kết quả",
    ready: "Sẵn sàng",
    running: "Đang chạy",
    completed: "Hoàn tất",
    failed: "Không đạt",
    cancelled: "Đã hủy",
    unknown: "Không rõ",
    jobs: "Tác vụ",
    systemStatus: "Trạng thái hệ thống",
    backend: "Backend",
    aiServer: "AI server",
    runRoot: "Thư mục chạy",
    fileRoots: "Vùng duyệt file",
    diskFree: "Ổ đĩa trống",
    runStorage: "Dung lượng kết quả",
    archive: "Tải gói ZIP",
    deleteJob: "Xóa tác vụ",
    validationRequired: "Vui lòng điền đủ thông tin bắt buộc.",
    invalidPrimer: "Primer chỉ được chứa A/T/G/C và mã IUPAC thoái hóa.",
    latest: "Lần chạy mới nhất",
    assays: "Bộ mồi",
    files: "Tệp",
    noJobs: "Chưa có lịch sử tác vụ.",
    selectJob: "Chọn tác vụ",
    routes: "Chức năng",
    localTitle: "Thiết kế từ tệp cục bộ",
    autoTitle: "Thiết kế tự động bằng từ khóa NCBI",
    aiTitle: "Thiết kế với AI",
    validateTitle: "Kiểm chứng bộ mồi đã biết",
    validationAction: "Chế độ thực hiện",
    validateMode: "Kiểm chứng mồi (Validate)",
    probeMode: "Thiết kế Probe cho mồi (Design Probes)",
    multiplexTitle: "Thiết kế PCR đa mồi",
    targetPath: "Bộ genome đích",
    backgroundPath: "Bộ genome nền",
    outputName: "Tên lần chạy",
    language: "Ngôn ngữ",
    runLocal: "Chạy thiết kế cục bộ",
    email: "Email NCBI",
    rememberEmail: "Ghi nhớ email NCBI",
    targetQueries: "Từ khóa đích",
    backgroundQueries: "Từ khóa nền",
    addTarget: "Thêm đích",
    addBackground: "Thêm nền",
    runAuto: "Chạy thiết kế tự động",
    query: "Truy vấn",
    size: "Kích thước tối thiểu (Mb)",
    count: "Số lượng",
    type: "Loại",
    remove: "Xóa",
    aiEndpoint: "Địa chỉ AI",
    aiModel: "Mô hình AI",
    aiAutoRun: "AI tự chạy khi đủ cấu hình",
    detectModels: "Dò model",
    modelServerOffline: "Chưa phát hiện model offline",
    modelProvider: "Nguồn model",
    assistantDock: "Trợ lý thường trực",
    assistantOpen: "Mở trợ lý",
    assistantClose: "Thu gọn",
    assistantPlaceholder: "Hỏi về app, PCR, qPCR, tối ưu wet-lab...",
    thinking: "AI đang suy nghĩ",
    browse: "Duyệt",
    choose: "Chọn",
    currentFolder: "Thư mục hiện tại",
    parentFolder: "Thư mục cha",
    fileBrowser: "Duyệt tệp và thư mục",
    folder: "Thư mục",
    file: "Tệp",
    openFolder: "Mở thư mục",
    chatPlaceholder: "VD: Thiết kế mồi cho Streptococcus suis, loại trừ các loài Streptococcus gần gũi",
    send: "Gửi",
    proposal: "Đề xuất của AI",
    runProposal: "Chạy đề xuất",
    blocked: "Cần bổ sung",
    params: "Tham số lõi",
    preflight: "Kiểm tra trước khi chạy",
    readyItem: "Sẵn sàng",
    missingItem: "Cần bổ sung",
    enabledItem: "Đang bật",
    disabledItem: "Đang tắt",
    primersCount: "Mồi",
    iupacWarning: "IUPAC tăng độ bao phủ chủng biến dị nhưng có thể tăng phản ứng chéo; luôn kiểm tra background trước wet-lab.",
    resultInterpretation: "Diễn giải nhanh kết quả",
    copy: "Copy",
    newChat: "Chat mới",
    copyAll: "Copy toàn bộ",
    copyPrimers: "Copy primer",
    bestCandidate: "Ứng viên tốt nhất",
    sensitivity: "Độ nhạy",
    specificity: "Độ đặc hiệu",
    amplicon: "Amplicon",
    askLatestResult: "Đánh giá kết quả mới nhất",
    askWetlab: "Tối ưu wet-lab",
    askDegenerate: "Thiết kế IUPAC",
    askMultiplex: "Tạo multiplex",
    restore: "Mặc định",
    stop: "Dừng",
    refresh: "Làm mới",
    logWaiting: "Đang chờ nhật ký pipeline...",
    noFinal: "Chưa có FINAL_ASSAY.csv.",
    noValidation: "Chưa có PCR_Advanced_Report.csv.",
    noMultiplex: "Chưa có MULTIPLEX_KITS.csv.",
    sourceMode: "Nguồn dữ liệu",
    analyze: "Phân tích",
    primerName: "Tên marker",
    fwd: "Mồi xuôi",
    rev: "Mồi ngược",
    probe: "Đoạn dò (Probe)",
    addPrimer: "Thêm mồi",
    maxMismatch: "Sai khác tối đa",
    maxLen: "Amplicon tối đa",
    extractSeq: "Xuất amplicon",
    runValidation: "Chạy kiểm chứng mồi",
    assayType: "Loại xét nghiệm",
    localTargets: "Thư mục đích cục bộ",
    addFolder: "Thêm thư mục",
    runMultiplex: "Chạy thiết kế đa mồi",
    existingFolders: "Thư mục kết quả",
    analyzeOnly: "Phân tích kết quả có sẵn",
    crossReactive: "Phản ứng chéo",
    jobLabel: "Tác vụ",
    sourceLabel: "Nguồn",
    outputLabel: "Đầu ra",
    userLabel: "Bạn",
    aiLabel: "AI",
    genome: "Genome",
    gene: "Gene",
    conventional: "PCR thường",
    unknownError: "Lỗi không xác định",
    min_sensitivity: "Độ nhạy tối thiểu",
    design_min_conservation: "Bảo tồn tối thiểu",
    design_max_candidates: "Số ứng viên tối đa",
    product_size_min: "Sản phẩm tối thiểu",
    product_size_max: "Sản phẩm tối đa",
    max_mismatch: "Sai khác tối đa",
    cpu_cores: "Số nhân CPU",
    enable_blast: "Bật BLAST",
    auto_relax_constraints: "Tự nới điều kiện",
    degenerate_primers: "Thiết kế mồi thoái hóa IUPAC",
    designSet: "Bộ mẫu thiết kế (Design Set)",
    validationSet: "Bộ mẫu kiểm chứng (Validation Set)",
    general: "Tham số chung",
    design_target_sampling_size: "Số lượng mẫu target (Design)",
    design_background_sampling_size: "Số lượng mẫu nền (Design)",
    validation_target_sampling_size: "Số lượng mẫu target (Validation)",
    validation_background_sampling_size: "Số lượng mẫu nền (Validation)",
    max_iupac_per_primer: "Số vị trí IUPAC tối đa/mồi",
    primer_length_min: "Chiều dài mồi tối thiểu",
    primer_length_max: "Chiều dài mồi tối đa",
    primer_tm_min: "Nhiệt độ nóng chảy tối thiểu",
    primer_tm_max: "Nhiệt độ nóng chảy tối đa",
    source_local: "Cục bộ",
    source_auto: "Tự động",
    source_upload: "Tải lên",
    source_validate_local: "Kiểm chứng cục bộ",
    source_validate_auto: "Kiểm chứng tự động",
    source_multiplex_local: "Đa mồi cục bộ",
    source_multiplex_auto: "Đa mồi tự động",
    source_multiplex_analyze: "Phân tích đa mồi",
    source_ai_local: "AI cục bộ",
    source_ai_auto: "AI tự động",
  },
  en: {
    app: "Rational Primer Design",
    subtitle: "qPCR assay design workspace",
    sharedConfig: "Shared configuration",
    aiConfig: "AI Expert",
    dashboard: "Dashboard",
    dashboardTitle: "Dashboard",
    local: "Local file",
    auto: "Auto keywords",
    ai: "AI chat",
    validate: "Validate primers",
    multiplex: "Multiplex",
    about: "About",
    aboutTitle: "About Rational Primer Design",
    aboutIntro:
      "Rational Primer Design is an integrated PCR/qPCR assay design, validation, and analysis workspace that combines an offline bioinformatics pipeline with a molecular diagnostics AI Expert.",
    aboutDesignTitle: "Controlled Assay Design",
    aboutDesignText:
      "The system supports design from local FASTA files, automatic NCBI keyword-based data retrieval, known-primer validation, and multiplex panel construction. Core parameters live in the sidebar so every task uses a consistent configuration.",
    aboutAiTitle: "PCR/qPCR AI Expert",
    aboutAiText:
      "The AI does more than generate runnable configurations. It can advise on target/background selection, primer/probe evaluation, in-silico result interpretation, dimer/hairpin/cross-reactivity risks, and wet-lab optimization plans.",
    aboutLocalTitle: "Local Data and Offline Models",
    aboutLocalText:
      "The package can work with local datasets, browse folders directly from the interface, auto-detect OpenAI-compatible/Ollama offline models, and preserve job history for result tracking.",
    aboutWorkflowTitle: "Primary Workflow",
    aboutWorkflowDesign: "Design singleplex or multiplex assays from target/background genomes.",
    aboutWorkflowValidate: "Validate known primers with advanced virtual PCR simulation.",
    aboutWorkflowReview: "Review logs, result tables, AI reports, and output files from the dashboard.",
    history: "History",
    monitor: "Monitor",
    results: "Results",
    ready: "Ready",
    running: "Running",
    completed: "Completed",
    failed: "Failed",
    cancelled: "Cancelled",
    unknown: "Unknown",
    jobs: "Jobs",
    systemStatus: "System status",
    backend: "Backend",
    aiServer: "AI server",
    runRoot: "Run folder",
    fileRoots: "File browser roots",
    diskFree: "Free disk",
    runStorage: "Result storage",
    archive: "Download ZIP",
    deleteJob: "Delete job",
    validationRequired: "Please fill in all required information.",
    invalidPrimer: "Primers may only contain A/T/G/C and IUPAC degenerate codes.",
    latest: "Latest run",
    assays: "Assays",
    files: "Files",
    noJobs: "No job history yet.",
    selectJob: "Select job",
    routes: "Dashboard",
    localTitle: "Design with local files",
    autoTitle: "Automatic design from NCBI keywords",
    aiTitle: "Design with AI chat",
    validateTitle: "Validate known primers",
    validationAction: "Action Mode",
    validateMode: "Validate Primers",
    probeMode: "Design Probes for Primers",
    multiplexTitle: "Multiplex PCR design",
    targetPath: "Target genomes",
    backgroundPath: "Background genomes",
    outputName: "Run name",
    language: "Language",
    runLocal: "Run local design",
    email: "NCBI email",
    rememberEmail: "Remember NCBI email",
    targetQueries: "Target keywords",
    backgroundQueries: "Background keywords",
    addTarget: "Add target",
    addBackground: "Add background",
    runAuto: "Run auto design",
    query: "Query",
    size: "Min size (Mb)",
    count: "Count",
    type: "Type",
    remove: "Remove",
    aiEndpoint: "AI endpoint",
    aiModel: "AI model",
    aiAutoRun: "AI auto-runs when configuration is complete",
    detectModels: "Detect models",
    modelServerOffline: "No offline model detected",
    modelProvider: "Model source",
    assistantDock: "Resident assistant",
    assistantOpen: "Open assistant",
    assistantClose: "Collapse",
    assistantPlaceholder: "Ask about the app, PCR, qPCR, wet-lab optimization...",
    thinking: "AI is thinking",
    browse: "Browse",
    choose: "Choose",
    currentFolder: "Current folder",
    parentFolder: "Parent folder",
    fileBrowser: "File and folder browser",
    folder: "Folder",
    file: "File",
    openFolder: "Open folder",
    chatPlaceholder: "Example: Design primers for Streptococcus suis and exclude close Streptococcus relatives",
    send: "Send",
    proposal: "AI proposal",
    runProposal: "Run proposal",
    blocked: "Needs input",
    params: "Core parameters",
    preflight: "Pre-flight check",
    readyItem: "Ready",
    missingItem: "Needs input",
    enabledItem: "Enabled",
    disabledItem: "Disabled",
    primersCount: "Primers",
    iupacWarning: "IUPAC increases coverage of variant strains but can increase cross-reactivity; always screen backgrounds before wet-lab use.",
    resultInterpretation: "Quick result interpretation",
    copy: "Copy",
    newChat: "New chat",
    copyAll: "Copy all",
    copyPrimers: "Copy primers",
    bestCandidate: "Best candidate",
    sensitivity: "Sensitivity",
    specificity: "Specificity",
    amplicon: "Amplicon",
    askLatestResult: "Review latest result",
    askWetlab: "Optimize wet-lab",
    askDegenerate: "Design IUPAC",
    askMultiplex: "Create multiplex",
    restore: "Defaults",
    stop: "Stop",
    refresh: "Refresh",
    logWaiting: "Waiting for pipeline logs...",
    noFinal: "No FINAL_ASSAY.csv yet.",
    noValidation: "No PCR_Advanced_Report.csv yet.",
    noMultiplex: "No MULTIPLEX_KITS.csv yet.",
    sourceMode: "Data source",
    analyze: "Analyze",
    primerName: "Marker name",
    fwd: "Forward",
    rev: "Reverse",
    probe: "Probe (Optional)",
    addPrimer: "Add primer",
    maxMismatch: "Max mismatch",
    maxLen: "Max amplicon",
    extractSeq: "Extract amplicon",
    runValidation: "Run primer validation",
    assayType: "Assay type",
    localTargets: "Local target folders",
    addFolder: "Add folder",
    runMultiplex: "Run multiplex",
    existingFolders: "Result folders",
    analyzeOnly: "Analyze existing",
    crossReactive: "Cross-reactive",
    jobLabel: "Job",
    sourceLabel: "Source",
    outputLabel: "Output",
    userLabel: "You",
    aiLabel: "AI",
    genome: "Genome",
    gene: "Gene",
    conventional: "Conventional",
    unknownError: "Unknown error",
    min_sensitivity: "Minimum sensitivity",
    design_min_conservation: "Minimum conservation",
    design_max_candidates: "Maximum candidates",
    product_size_min: "Minimum product size",
    product_size_max: "Maximum product size",
    max_mismatch: "Maximum mismatch",
    cpu_cores: "CPU cores",
    enable_blast: "Enable BLAST",
    auto_relax_constraints: "Auto-relax constraints",
    degenerate_primers: "Design IUPAC degenerate primers",
    designSet: "Design Set",
    validationSet: "Validation Set",
    general: "General",
    design_target_sampling_size: "Design target sample size",
    design_background_sampling_size: "Design background sample size",
    validation_target_sampling_size: "Validation target sample size",
    validation_background_sampling_size: "Validation background sample size",
    max_iupac_per_primer: "Max IUPAC positions per primer",
    primer_length_min: "Min primer length",
    primer_length_max: "Max primer length",
    primer_tm_min: "Min primer Tm",
    primer_tm_max: "Max primer Tm",
    source_local: "Local",
    source_auto: "Automatic",
    source_upload: "Upload",
    source_validate_local: "Local validation",
    source_validate_auto: "Automatic validation",
    source_multiplex_local: "Local multiplex",
    source_multiplex_auto: "Automatic multiplex",
    source_multiplex_analyze: "Multiplex analysis",
    source_ai_local: "AI local",
    source_ai_auto: "AI automatic",
  },
};

const PARAM_GROUPS = ["designSet", "validationSet", "general"] as const;

type ParamField = {
  key: string;
  type: "number" | "boolean";
  group: (typeof PARAM_GROUPS)[number];
  min?: number;
  max?: number;
  step?: number;
};

const PARAM_FIELDS: ParamField[] = [
  { key: "design_target_sampling_size", type: "number", group: "designSet", min: 0, max: 500, step: 1 },
  { key: "design_background_sampling_size", type: "number", group: "designSet", min: 0, max: 500, step: 1 },
  { key: "primer_length_min", type: "number", group: "designSet", min: 10, max: 50, step: 1 },
  { key: "primer_length_max", type: "number", group: "designSet", min: 10, max: 50, step: 1 },
  { key: "primer_tm_min", type: "number", group: "designSet", min: 40, max: 80, step: 0.5 },
  { key: "primer_tm_max", type: "number", group: "designSet", min: 40, max: 80, step: 0.5 },
  { key: "design_min_conservation", type: "number", group: "designSet", min: 0.5, max: 1, step: 0.01 },
  { key: "design_max_candidates", type: "number", group: "designSet", min: 1, max: 500, step: 1 },
  { key: "validation_target_sampling_size", type: "number", group: "validationSet", min: 0, max: 500, step: 1 },
  { key: "validation_background_sampling_size", type: "number", group: "validationSet", min: 0, max: 500, step: 1 },
  { key: "min_sensitivity", type: "number", group: "validationSet", min: 50, max: 100, step: 0.1 },
  { key: "max_mismatch", type: "number", group: "validationSet", min: 0, max: 8, step: 1 },
  { key: "product_size_min", type: "number", group: "general", min: 50, max: 1000, step: 1 },
  { key: "product_size_max", type: "number", group: "general", min: 50, max: 2000, step: 1 },
  { key: "cpu_cores", type: "number", group: "general", min: 0, max: 128, step: 1 },
  { key: "enable_blast", type: "boolean", group: "general" },
  { key: "auto_relax_constraints", type: "boolean", group: "general" },
  { key: "degenerate_primers", type: "boolean", group: "general" },
  { key: "max_iupac_per_primer", type: "number", group: "general", min: 0, max: 10, step: 1 },
];

const NAV: Array<{ view: WorkspaceView; href: string; icon: string }> = [
  { view: "dashboard", href: "/dashboard", icon: "▤" },
  { view: "local", href: "/design/local", icon: "⌂" },
  { view: "auto", href: "/design/auto", icon: "◎" },
  { view: "ai", href: "/ai", icon: "◆" },
  { view: "validate", href: "/validate", icon: "✓" },
  { view: "multiplex", href: "/multiplex", icon: "⊞" },
  { view: "about", href: "/about", icon: "ⓘ" },
];

export function Dashboard({ view = "dashboard" }: { view?: WorkspaceView }) {
  const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;
  const [lang, setLang] = useState<Lang>("en");
  const [parameters, setParameters] = useState<Parameters>({});
  const [jobs, setJobs] = useState<Job[]>([]);
  const [job, setJob] = useState<Job | null>(null);
  const [logs, setLogs] = useState("");
  const [results, setResults] = useState<ResultPayload | null>(null);
  const [status, setStatus] = useState<StatusPayload | null>(null);
  const [error, setError] = useState("");
  const [isBusy, setIsBusy] = useState(false);

  const [localTarget, setLocalTarget] = useState("");
  const [localBackground, setLocalBackground] = useState("");
  const [outputName, setOutputName] = useState("primer-design-run");
  const [localMode, setLocalMode] = useState<"path" | "upload">("path");
  const [uploadTargetFiles, setUploadTargetFiles] = useState<File[]>([]);
  const [uploadBackgroundFiles, setUploadBackgroundFiles] = useState<File[]>([]);

  const [email, setEmail] = useState("");
  const [rememberEmail, setRememberEmail] = useState(false);
  const [targets, setTargets] = useState<QueryItem[]>([{ query: "", size: 0, count: 500, type: "genome" }]);
  const [backgrounds, setBackgrounds] = useState<QueryItem[]>([{ query: "", size: 0, count: 50, type: "genome" }]);

  const [validationMode, setValidationMode] = useState<DataMode>("local");
  const [validationAction, setValidationAction] = useState<"validate" | "probe">("validate");
  const [validationPrimers, setValidationPrimers] = useState<PrimerPair[]>([{ name: "M1", fwd: "", rev: "", probe: "" }]);
  const [validationTarget, setValidationTarget] = useState("test_data/target");
  const [validationBackground, setValidationBackground] = useState("test_data/background");
  const [validationTargets, setValidationTargets] = useState<QueryItem[]>([{ query: "", size: 0, count: 50, type: "genome" }]);
  const [validationBackgrounds, setValidationBackgrounds] = useState<QueryItem[]>([{ query: "", size: 0, count: 10, type: "genome" }]);
  const [maxMismatch, setMaxMismatch] = useState(4);
  const [maxLen, setMaxLen] = useState(1500);
  const [extractSeq, setExtractSeq] = useState(false);

  const [multiplexMode, setMultiplexMode] = useState<MultiplexMode>("local");
  const [assayType, setAssayType] = useState<"qPCR" | "Conventional">("qPCR");
  const [localMultiplexTargets, setLocalMultiplexTargets] = useState<string[]>(["test_data/target", ""]);
  const [multiplexBackground, setMultiplexBackground] = useState("test_data/background");
  const [multiplexTargets, setMultiplexTargets] = useState<QueryItem[]>([
    { query: "", size: 0, count: 500, type: "genome" },
    { query: "", size: 0, count: 500, type: "genome" },
  ]);
  const [multiplexBackgrounds, setMultiplexBackgrounds] = useState<QueryItem[]>([{ query: "", size: 0, count: 50, type: "genome" }]);
  const [multiplexFolders, setMultiplexFolders] = useState<string[]>(["", ""]);

  const [aiBaseUrl, setAiBaseUrl] = useState("http://localhost:11434/v1");
  const [aiModel, setAiModel] = useState("llama3");
  const [offlineModels, setOfflineModels] = useState<OfflineModelsResponse | null>(null);
  const [aiAutoRun, setAiAutoRun] = useState(true);
  const [chatInput, setChatInput] = useState("");
  const [chat, setChat] = useState<ChatMessage[]>([
    {
      role: "assistant",
      content:
        "Tôi là AI Expert. Hãy nói target, background, bộ primer hoặc yêu cầu multiplex; tôi sẽ đề xuất cấu hình và có thể tự chạy pipeline.",
    },
  ]);
  const [proposal, setProposal] = useState<Record<string, unknown> | null>(null);
  const [blockedReason, setBlockedReason] = useState("");
  const [assistantOpen, setAssistantOpen] = useState(false);
  const [assistantInput, setAssistantInput] = useState("");
  const [assistantBusy, setAssistantBusy] = useState(false);
  const [assistantChat, setAssistantChat] = useState<ChatMessage[]>([
    {
      role: "assistant",
      content:
        "Tôi luôn ở đây để trả lời về app, thiết kế primer/PCR/qPCR, đánh giá kết quả và tối ưu wet-lab. Tôi sẽ không tự chạy pipeline từ khung chat nhỏ này.",
    },
  ]);

  function resetAssistantChat() {
    setAssistantChat([
      {
        role: "assistant",
        content:
          "Tôi luôn ở đây để trả lời về app, thiết kế primer/PCR/qPCR, đánh giá kết quả và tối ưu wet-lab. Tôi sẽ không tự chạy pipeline từ khung chat nhỏ này.",
      },
    ]);
    setProposal(null);
    setBlockedReason("");
  }

  function resetChat() {
    setChat([
      {
        role: "assistant",
        content:
          "Tôi là AI Expert. Hãy nói target, background, bộ primer hoặc yêu cầu multiplex; tôi sẽ đề xuất cấu hình và có thể tự chạy pipeline.",
      },
    ]);
    setProposal(null);
    setBlockedReason("");
  }
  const [fileBrowser, setFileBrowser] = useState<FileBrowserState | null>(null);

  const tx = I18N[lang];
  const title = viewTitle(view, tx);
  const statusText = job ? tx[job.status] ?? job.status : tx.ready;
  const crossSummary = results?.summary?.cross_contamination_summary as { total_cross_reactive_strains?: number } | undefined;
  const shellClassName = view === "ai" ? "site-shell ai-screen" : "site-shell";

  useEffect(() => {
    void loadDefaults();
    void refreshJobs();
    const savedLanguage = window.localStorage.getItem("rpd-language");
    if (savedLanguage === "vi" || savedLanguage === "en") setLang(savedLanguage);
    const savedRememberEmail = window.localStorage.getItem("rpd-remember-email") === "1";
    setRememberEmail(savedRememberEmail);
    if (savedRememberEmail) setEmail(window.localStorage.getItem("rpd-ncbi-email") ?? "");
    const savedAiBaseUrl = window.localStorage.getItem("rpd-ai-base-url") ?? "http://localhost:11434/v1";
    const savedAiModel = window.localStorage.getItem("rpd-ai-model") ?? "llama3";
    setAiBaseUrl(savedAiBaseUrl);
    setAiModel(savedAiModel);
    setOutputName(window.localStorage.getItem("rpd-output-name") ?? "primer-design-run");
    setLocalTarget(window.localStorage.getItem("rpd-local-target") ?? "");
    setLocalBackground(window.localStorage.getItem("rpd-local-background") ?? "");
    setValidationTarget(window.localStorage.getItem("rpd-validation-target") ?? "");
    setValidationBackground(window.localStorage.getItem("rpd-validation-background") ?? "");
    setMultiplexBackground(window.localStorage.getItem("rpd-multiplex-background") ?? "");
    void detectOfflineModels(savedAiBaseUrl, savedAiModel);
  }, []);

  useEffect(() => {
    window.localStorage.setItem("rpd-remember-email", rememberEmail ? "1" : "0");
    if (rememberEmail) {
      window.localStorage.setItem("rpd-ncbi-email", email);
    } else {
      window.localStorage.removeItem("rpd-ncbi-email");
    }
  }, [rememberEmail, email]);

  useEffect(() => {
    setChat((current) => {
      if (current.length !== 1 || current[0].role !== "assistant") return current;
      return [
        {
          role: "assistant",
          content:
            lang === "vi"
              ? "Tôi là AI Expert. Hãy nói target, background, bộ primer hoặc yêu cầu multiplex; tôi sẽ đề xuất cấu hình và có thể tự chạy pipeline."
              : "I am the AI Expert. Tell me the target, background, primer set, or multiplex request; I will propose a configuration and can run the pipeline.",
        },
      ];
    });
  }, [lang]);

  useEffect(() => {
    setAssistantChat((current) => {
      if (current.length !== 1 || current[0].role !== "assistant") return current;
      return [
        {
          role: "assistant",
          content:
            lang === "vi"
              ? "Tôi luôn ở đây để trả lời về app, thiết kế primer/PCR/qPCR, đánh giá kết quả và tối ưu wet-lab. Tôi sẽ không tự chạy pipeline từ khung chat nhỏ này."
              : "I am always here for questions about the app, primer/PCR/qPCR design, result interpretation, and wet-lab optimization. I will not auto-run pipelines from this small chat.",
        },
      ];
    });
  }, [lang]);

  const didMount = useRef(false);
  useEffect(() => {
    if (!didMount.current) { didMount.current = true; return; }
    const timer = window.setTimeout(() => {
      void detectOfflineModels();
      void refreshStatus();
    }, 400);
    return () => window.clearTimeout(timer);
  }, [aiBaseUrl]);

  useEffect(() => {
    window.localStorage.setItem("rpd-ai-base-url", aiBaseUrl);
    window.localStorage.setItem("rpd-ai-model", aiModel);
    window.localStorage.setItem("rpd-output-name", outputName);
    window.localStorage.setItem("rpd-local-target", localTarget);
    window.localStorage.setItem("rpd-local-background", localBackground);
    window.localStorage.setItem("rpd-validation-target", validationTarget);
    window.localStorage.setItem("rpd-validation-background", validationBackground);
    window.localStorage.setItem("rpd-multiplex-background", multiplexBackground);
  }, [aiBaseUrl, aiModel, outputName, localTarget, localBackground, validationTarget, validationBackground, multiplexBackground]);

  useEffect(() => {
    if (!job) return;
    void refreshJob(job.id);
    const timer = window.setInterval(() => void refreshJob(job.id), job.status === "running" ? 2000 : 6000);
    return () => window.clearInterval(timer);
  }, [job?.id, job?.status]);

  async function api<T>(path: string, init?: RequestInit): Promise<T> {
    const response = await fetch(`${API_BASE}${path}`, init);
    if (!response.ok) {
      let detail: unknown = response.statusText;
      try {
        const payload = (await response.json()) as Record<string, unknown>;
        detail = payload.detail ?? detail;
      } catch {
        detail = await response.text();
      }
      throw new Error(Array.isArray(detail) ? detail.map((d: unknown) => typeof d === "object" && d ? String((d as Record<string, unknown>).msg ?? JSON.stringify(d)) : String(d)).join("; ") : String(detail));
    }
    return response.json() as Promise<T>;
  }

  async function loadDefaults() {
    try {
      setParameters(await api<Parameters>("/api/default-parameters"));
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    }
  }

  async function refreshJobs() {
    try {
      const data = await api<Job[]>("/api/jobs");
      setJobs(data);
      const history = await api<HistoryPayload>("/api/history");
    } catch {
      setJobs([]);
    }
  }

  async function refreshJob(jobId: string) {
    try {
      const nextJob = await api<Job>(`/api/jobs/${jobId}`);
      setJob(nextJob);
      const logPayload = await api<{ logs: string }>(`/api/jobs/${jobId}/logs`);
      setLogs(logPayload.logs);
      setResults(await api<ResultPayload>(`/api/jobs/${jobId}/results`));
      void refreshJobs();
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    }
  }

  async function refreshStatus(baseUrlOverride?: string) {
    try {
      setStatus(await api<StatusPayload>(`/api/status?ai_base_url=${encodeURIComponent(baseUrlOverride ?? aiBaseUrl)}`));
    } catch {
      setStatus(null);
    }
  }

  async function detectOfflineModels(baseUrlOverride?: string, currentModel?: string) {
    const url = baseUrlOverride ?? aiBaseUrl;
    const model = currentModel ?? aiModel;
    try {
      const response = await api<OfflineModelsResponse>(`/api/ai/models?base_url=${encodeURIComponent(url)}`);
      setOfflineModels(response);
      if (response.models.length > 0 && (!model || model === "llama3")) {
        setAiModel(response.models[0]);
      }
      void refreshStatus(url);
    } catch (err) {
      setOfflineModels({
        base_url: url,
        provider: "unknown",
        available: false,
        models: [],
        error: messageOf(err, tx.unknownError),
      });
    }
  }

  function setParam(key: string, value: string | number | boolean) {
    setParameters((current) => ({ ...current, [key]: value }));
  }

  function changeLanguage(value: Lang) {
    setLang(value);
    window.localStorage.setItem("rpd-language", value);
  }

  function requestBrowse(path: string, kind: BrowseKind, onSelect: (path: string) => void) {
    setFileBrowser({ open: true, path: path || "~", kind, onSelect });
  }

  async function startLocal(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (!localTarget.trim() || !localBackground.trim()) return setError(tx.validationRequired);
    await createJob("/api/jobs/local", {
      target_path: localTarget,
      background_path: localBackground,
      output_name: outputName,
      parameters,
      language: lang,
      ai_base_url: aiBaseUrl,
      ai_model: aiModel,
    });
  }

  async function startUpload(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (uploadTargetFiles.length === 0 || uploadBackgroundFiles.length === 0) return setError(tx.validationRequired);
    setIsBusy(true);
    setLogs("");
    setResults(null);
    setError("");
    try {
      const formData = new FormData();
      for (const file of uploadTargetFiles) formData.append("target_files", file);
      for (const file of uploadBackgroundFiles) formData.append("background_files", file);
      if (outputName.trim()) formData.append("output_name", outputName.trim());
      formData.append("parameters", JSON.stringify(parameters));
      formData.append("language", lang);
      formData.append("ai_base_url", aiBaseUrl);
      formData.append("ai_model", aiModel);
      const created = await api<Job>("/api/jobs/upload", {
        method: "POST",
        body: formData,
      });
      setJob(created);
      setUploadTargetFiles([]);
      setUploadBackgroundFiles([]);
      void refreshJobs();
      void refreshStatus();
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    } finally {
      setIsBusy(false);
    }
  }

  async function startAuto(event?: FormEvent<HTMLFormElement>) {
    event?.preventDefault();
    if (!validEmail(email) || nonEmptyQueries(targets).length === 0 || nonEmptyQueries(backgrounds).length === 0) return setError(tx.validationRequired);
    await createJob("/api/jobs/auto", {
      email,
      targets: nonEmptyQueries(targets),
      backgrounds: nonEmptyQueries(backgrounds),
      output_name: outputName,
      parameters,
      language: lang,
      ai_base_url: aiBaseUrl,
      ai_model: aiModel,
    });
  }

  async function startValidation(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const cleanPrimers = validationPrimers.filter((item) => item.name.trim() && item.fwd.trim() && item.rev.trim());
    if (cleanPrimers.length === 0) return setError(tx.validationRequired);
    if (!cleanPrimers.every(validPrimerPair)) return setError(tx.invalidPrimer);
    if (validationMode === "local" && !validationTarget.trim()) return setError(tx.validationRequired);
    if (validationMode === "auto" && (!validEmail(email) || nonEmptyQueries(validationTargets).length === 0)) return setError(tx.validationRequired);
    const base = {
      primers: cleanPrimers,
      output_name: outputName,
      extract_sequence: extractSeq,
      max_mismatch: maxMismatch,
      max_len: maxLen,
      action: validationAction === "probe" ? "design-probes" : "validate",
    };
    await createJob(
      "/api/jobs/validate",
      validationMode === "local"
        ? { ...base, target_path: validationTarget, background_path: validationBackground }
        : {
            ...base,
            email,
            targets: nonEmptyQueries(validationTargets),
            backgrounds: nonEmptyQueries(validationBackgrounds),
          },
    );
  }

  async function startMultiplex(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (multiplexMode === "local") {
      if (localMultiplexTargets.filter((item) => item.trim()).length < 2 || !multiplexBackground.trim()) return setError(tx.validationRequired);
      await createJob("/api/jobs/multiplex/local", {
        target_paths: localMultiplexTargets.filter((item) => item.trim()),
        background_path: multiplexBackground,
        output_name: outputName,
        parameters,
        language: lang,
        assay_type: assayType,
        ai_base_url: aiBaseUrl,
        ai_model: aiModel,
      });
    } else if (multiplexMode === "auto") {
      if (!validEmail(email) || nonEmptyQueries(multiplexTargets).length < 2) return setError(tx.validationRequired);
      await createJob("/api/jobs/multiplex/auto", {
        email,
        targets: nonEmptyQueries(multiplexTargets),
        backgrounds: nonEmptyQueries(multiplexBackgrounds),
        output_name: outputName,
        parameters,
        language: lang,
        assay_type: assayType,
        ai_base_url: aiBaseUrl,
        ai_model: aiModel,
      });
    } else {
      if (multiplexFolders.filter((item) => item.trim()).length < 2) return setError(tx.validationRequired);
      await createJob("/api/jobs/multiplex/analyze", {
        folders: multiplexFolders.filter((item) => item.trim()),
        output_name: outputName,
        language: lang,
        assay_type: assayType,
        ai_base_url: aiBaseUrl,
        ai_model: aiModel,
      });
    }
  }

  async function createJob(path: string, body: Record<string, unknown>) {
    setError("");
    setIsBusy(true);
    setLogs("");
    setResults(null);
    try {
      const created = await api<Job>(path, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      setJob(created);
      void refreshJobs();
      void refreshStatus();
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    } finally {
      setIsBusy(false);
    }
  }

  async function sendChat(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (!chatInput.trim()) return;
    const nextMessages: ChatMessage[] = [...chat, { role: "user", content: chatInput.trim() }];
    setChat(nextMessages);
    setChatInput("");
    setBlockedReason("");
    setError("");
    setIsBusy(true);
    try {
      const response = await api<ChatResponse>("/api/ai/chat", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          messages: nextMessages,
          language: lang,
          ai_base_url: aiBaseUrl,
          ai_model: aiModel,
          email,
          output_name: outputName,
          parameters,
          auto_run: aiAutoRun,
          expert_context: results,
        }),
      });
      setChat((current) => [...current, { role: "assistant", content: response.reply || response.raw_reply }]);
      setProposal(response.proposal);
      setBlockedReason(response.blocked_reason ?? "");
      if (response.job) {
        setJob(response.job);
        void refreshJobs();
      }
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    } finally {
      setIsBusy(false);
    }
  }

  async function sendAssistantChat(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (!assistantInput.trim()) return;
    const nextMessages: ChatMessage[] = [
      ...assistantChat,
      {
        role: "user",
        content: assistantInput.trim(),
      },
    ];
    setAssistantChat(nextMessages);
    setAssistantInput("");
    setAssistantBusy(true);
    setError("");
    try {
      const response = await api<ChatResponse>("/api/ai/chat", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          messages: nextMessages,
          language: lang,
          ai_base_url: aiBaseUrl,
          ai_model: aiModel,
          email,
          output_name: outputName,
          parameters,
          auto_run: false,
          expert_context: {
            current_view: view,
            selected_job: job,
            current_results: results,
            app_capabilities: ["dashboard_history", "local_design", "auto_ncbi_design", "ai_design", "known_primer_validation", "multiplex_design"],
          },
        }),
      });
      setAssistantChat((current) => [...current, { role: "assistant", content: response.reply || response.raw_reply }]);
    } catch (err) {
      setAssistantChat((current) => [...current, { role: "assistant", content: messageOf(err, tx.unknownError) }]);
    } finally {
      setAssistantBusy(false);
    }
  }

  async function runProposal() {
    if (!proposal) return;
    const action = String(proposal.action ?? "");
    const needsEmail = action === "propose_design" || action === "propose_validation" || action === "propose_multiplex";
    if (needsEmail && !validEmail(email)) {
      return setError(tx.email + " " + tx.validationRequired);
    }
    if (action === "propose_local_design" && (!String(proposal.local_target ?? "").trim() || !String(proposal.local_background ?? "").trim())) {
      return setError(tx.validationRequired);
    }
    setProposal(null);
    if (action === "propose_local_design") {
      await createJob("/api/jobs/local", {
        target_path: String(proposal.local_target ?? "").trim(),
        background_path: String(proposal.local_background ?? "").trim(),
        output_name: outputName,
        parameters: proposalParameters(proposal, parameters),
        language: lang,
        ai_base_url: aiBaseUrl,
        ai_model: aiModel,
      });
    } else if (action === "propose_design") {
      await createJob("/api/jobs/auto", {
        email,
        targets: proposalList(proposal, "target", "targets"),
        backgrounds: proposalList(proposal, "background", "backgrounds"),
        output_name: outputName,
        parameters: proposalParameters(proposal, parameters),
        language: lang,
        ai_base_url: aiBaseUrl,
        ai_model: aiModel,
      });
    } else if (action === "propose_validation") {
      await createJob("/api/jobs/validate", {
        email,
        primers: proposalList(proposal, "primers"),
        targets: proposalList(proposal, "target", "targets"),
        backgrounds: proposalList(proposal, "background", "backgrounds"),
        output_name: outputName,
        max_mismatch: Number(proposal.max_mismatch ?? 4),
        max_len: Number(proposal.max_len ?? 1500),
      });
    } else if (action === "propose_multiplex") {
      await createJob("/api/jobs/multiplex/auto", {
        email,
        targets: proposalList(proposal, "targets", "target"),
        backgrounds: proposalList(proposal, "background", "backgrounds"),
        output_name: outputName,
        parameters: proposalParameters(proposal, parameters),
        language: lang,
        assay_type: String(proposal.assay_type ?? "qPCR"),
        ai_base_url: aiBaseUrl,
        ai_model: aiModel,
      });
    } else {
      setBlockedReason(`${tx.blocked}: ${action}`);
    }
  }

  async function cancelCurrentJob() {
    if (!job) return;
    try {
      setJob(await api<Job>(`/api/jobs/${job.id}/cancel`, { method: "POST" }));
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    }
  }

  async function deleteCurrentJob() {
    if (!job) return;
    try {
      await api<{ deleted: boolean }>(`/api/jobs/${job.id}`, { method: "DELETE" });
      setJob(null);
      setLogs("");
      setResults(null);
      void refreshJobs();
      void refreshStatus();
    } catch (err) {
      setError(messageOf(err, tx.unknownError));
    }
  }

  const terminal = useMemo(() => logs || (job ? tx.logWaiting : tx.ready), [job, logs, tx.logWaiting, tx.ready]);

  return (
    <main className={shellClassName}>
      <header className="site-header">
        <div className="brand">
          <div className="brand-mark">RP</div>
          <div>
            <h1>{tx.app}</h1>
            <p>{tx.subtitle}</p>
          </div>
        </div>
        <nav className="top-nav" aria-label={tx.routes}>
          {NAV.filter((n) => n.view !== "about").map((item) => (
            <Link key={item.view} href={item.href} className={view === item.view ? "top-nav-link active" : "top-nav-link"}>
              <span className="nav-icon" aria-hidden="true">{item.icon}</span>
              <span className="nav-label">{tx[item.view]}</span>
            </Link>
          ))}
          <Link href="/about" className={`top-nav-link about-link ${view === "about" ? "active" : ""}`}>
            <span className="nav-icon" aria-hidden="true">ⓘ</span>
          </Link>
        </nav>
        <div className="status-strip">
          <span className={`dot ${job?.status ?? "idle"}`} />
          <span>{statusText}</span>
        </div>
      </header>

      <div className={"site-body" + (view === "dashboard" ? " site-body--full" : "")}>
        {view !== "dashboard" ? (
          <aside className="config-sidebar" aria-label={tx.sharedConfig}>
            <SharedConfigPanel
              tx={tx}
              lang={lang}
              changeLanguage={changeLanguage}
              outputName={outputName}
              setOutputName={setOutputName}
              email={email}
              setEmail={setEmail}
              rememberEmail={rememberEmail}
              setRememberEmail={setRememberEmail}
              aiBaseUrl={aiBaseUrl}
              setAiBaseUrl={setAiBaseUrl}
              aiModel={aiModel}
              setAiModel={setAiModel}
              offlineModels={offlineModels}
              detectOfflineModels={detectOfflineModels}
              aiAutoRun={aiAutoRun}
              setAiAutoRun={setAiAutoRun}
              assayType={assayType}
              setAssayType={setAssayType}
              parameters={parameters}
              setParam={setParam}
              restore={() => void loadDefaults()}
            />
            <SystemStatusPanel tx={tx} status={status} refresh={() => void refreshStatus()} />
            <RecentJobs tx={tx} jobs={jobs} select={(item) => setJob(item)} />
          </aside>
        ) : null}

        <section className="workspace">
          <header className="topbar">
            <div>
              {view !== "dashboard" ? <div className="eyebrow">{tx.routes}</div> : null}
              <h2>{title}</h2>
            </div>
          </header>

          {error ? <div className="alert">{error}</div> : null}

          {view !== "about" ? (
            <section className="dashboard-grid">
              {view === "dashboard" ? (
                <>
                  <Metric label={tx.jobs} value={String(jobs.length)} icon="▤" />
                  <Metric label={tx.latest} value={job ? job.id.slice(-12) : "-"} icon="◷" />
                  <Metric label={tx.assays} value={String(results?.final_assays.row_count ?? results?.multiplex_kits.row_count ?? 0)} icon="⊞" />
                  <Metric label={tx.crossReactive} value={String(crossSummary?.total_cross_reactive_strains ?? 0)} icon="⚠" />
                </>
              ) : view === "local" || view === "auto" ? (
                <>
                  <Metric label={tx.jobs} value={String(jobs.length)} icon="▤" />
                  <Metric label={tx.outputName} value={outputName || "-"} icon="✎" />
                  <Metric label={tx.targetPath} value={view === "local" ? localTarget : (targets[0]?.query || "-")} icon="⌂" />
                  <Metric label={tx.assays} value={String(results?.final_assays.row_count ?? 0)} icon="⊞" />
                </>
              ) : view === "validate" ? (
                <>
                  <Metric label={tx.jobs} value={String(jobs.length)} icon="▤" />
                  <Metric label={tx.primersCount || "Primers"} value={String(validationPrimers.filter((p) => p.fwd.trim()).length)} icon="✓" />
                  <Metric label={tx.maxMismatch} value={String(maxMismatch)} icon="◷" />
                  <Metric label={tx.assays} value={String(results?.known_primer_validation.row_count ?? 0)} icon="⊞" />
                </>
              ) : view === "multiplex" ? (
                <>
                  <Metric label={tx.jobs} value={String(jobs.length)} icon="▤" />
                  <Metric label={tx.assayType} value={assayType} icon="⊞" />
                  <Metric label={tx.sourceMode} value={tx[`source_multiplex_${multiplexMode}`] || multiplexMode} icon="◎" />
                  <Metric label={tx.assays} value={String(results?.multiplex_kits.row_count ?? 0)} icon="⊞" />
                </>
              ) : view === "ai" ? (
                <>
                  <Metric label={tx.jobs} value={String(jobs.length)} icon="▤" />
                  <Metric label={tx.aiModel} value={aiModel || "-"} icon="◆" />
                  <Metric label={tx.aiEndpoint || "Endpoint"} value={aiBaseUrl ? aiBaseUrl.replace(/^https?:\/\//, "") : "-"} icon="◎" />
                  <Metric label={tx.assays} value={String(results?.final_assays.row_count ?? 0)} icon="⊞" />
                </>
              ) : null}
            </section>
          ) : null}

          {view === "about" ? (
            <AboutPage tx={tx} />
          ) : view === "dashboard" ? (
            <section className="main-grid">
              <div className="stack">
                <DashboardHome tx={tx} jobs={jobs} selectedJob={job} select={(item) => setJob(item)} />
                <ResultsPanel tx={tx} results={results} requestBrowse={requestBrowse} />
              </div>
              <MonitorPanel tx={tx} job={job} terminal={terminal} cancel={cancelCurrentJob} deleteJob={deleteCurrentJob} refresh={() => job && void refreshJob(job.id)} />
            </section>
          ) : (
            <section className={view === "ai" ? "main-grid ai-main-grid" : "main-grid"}>
              <div className="stack">
                {view === "local" ? (
                  <LocalMode
                    tx={tx}
                    isElectron={isElectron}
                    target={localTarget}
                    background={localBackground}
                    setTarget={setLocalTarget}
                    setBackground={setLocalBackground}
                    outputName={outputName}
                    setOutputName={setOutputName}
                    requestBrowse={requestBrowse}
                    parameters={parameters}
                    onSubmit={startLocal}
                    onUploadSubmit={startUpload}
                    isBusy={isBusy}
                    uploadTargetFiles={uploadTargetFiles}
                    uploadBackgroundFiles={uploadBackgroundFiles}
                    setUploadTargetFiles={setUploadTargetFiles}
                    setUploadBackgroundFiles={setUploadBackgroundFiles}
                  />
                ) : null}
                {view === "auto" ? (
                  <AutoMode
                    tx={tx}
                    email={email}
                    setEmail={setEmail}
                    targets={targets}
                    backgrounds={backgrounds}
                    setTargets={setTargets}
                    setBackgrounds={setBackgrounds}
                    outputName={outputName}
                    setOutputName={setOutputName}
                    parameters={parameters}
                    onSubmit={startAuto}
                    isBusy={isBusy}
                  />
                ) : null}
                {view === "ai" ? (
                    <AiMode
                      tx={tx}
                      email={email}
                      setEmail={setEmail}
                      aiBaseUrl={aiBaseUrl}
                      setAiBaseUrl={setAiBaseUrl}
                      aiModel={aiModel}
                      setAiModel={setAiModel}
                      aiAutoRun={aiAutoRun}
                      setAiAutoRun={setAiAutoRun}
                      lang={lang}
                      chat={chat}
                      chatInput={chatInput}
                      setChatInput={setChatInput}
                      sendChat={sendChat}
                      proposal={proposal}
                      blockedReason={blockedReason}
                      runProposal={runProposal}
                      isBusy={isBusy}
                      onNewChat={resetChat}
                    />
                ) : null}
                {view === "validate" ? (
                  <ValidationMode
                    tx={tx}
                    mode={validationMode}
                    setMode={setValidationMode}
                    action={validationAction}
                    setAction={setValidationAction}
                    primers={validationPrimers}
                    setPrimers={setValidationPrimers}
                    targetPath={validationTarget}
                    backgroundPath={validationBackground}
                    setTargetPath={setValidationTarget}
                    setBackgroundPath={setValidationBackground}
                    requestBrowse={requestBrowse}
                    email={email}
                    setEmail={setEmail}
                    targets={validationTargets}
                    backgrounds={validationBackgrounds}
                    setTargets={setValidationTargets}
                    setBackgrounds={setValidationBackgrounds}
                    outputName={outputName}
                    setOutputName={setOutputName}
                    maxMismatch={maxMismatch}
                    setMaxMismatch={setMaxMismatch}
                    maxLen={maxLen}
                    setMaxLen={setMaxLen}
                    extractSeq={extractSeq}
                    setExtractSeq={setExtractSeq}
                    parameters={parameters}
                    onSubmit={startValidation}
                    isBusy={isBusy}
                  />
                ) : null}
                {view === "multiplex" ? (
                  <MultiplexModePanel
                    tx={tx}
                    mode={multiplexMode}
                    setMode={setMultiplexMode}
                    assayType={assayType}
                    setAssayType={setAssayType}
                    localTargets={localMultiplexTargets}
                    setLocalTargets={setLocalMultiplexTargets}
                    backgroundPath={multiplexBackground}
                    setBackgroundPath={setMultiplexBackground}
                    requestBrowse={requestBrowse}
                    email={email}
                    setEmail={setEmail}
                    targets={multiplexTargets}
                    backgrounds={multiplexBackgrounds}
                    setTargets={setMultiplexTargets}
                    setBackgrounds={setMultiplexBackgrounds}
                    folders={multiplexFolders}
                    setFolders={setMultiplexFolders}
                    outputName={outputName}
                    setOutputName={setOutputName}
                    parameters={parameters}
                    onSubmit={startMultiplex}
                    isBusy={isBusy}
                  />
                ) : null}
              </div>

              <div className="stack">
                <MonitorPanel tx={tx} job={job} terminal={terminal} cancel={cancelCurrentJob} deleteJob={deleteCurrentJob} refresh={() => job && void refreshJob(job.id)} />
                <ResultsPanel tx={tx} results={results} requestBrowse={requestBrowse} />
              </div>
            </section>
          )}
        </section>
      </div>
      <FloatingAssistant
        tx={tx}
        open={assistantOpen}
        setOpen={setAssistantOpen}
        messages={assistantChat}
        input={assistantInput}
        setInput={setAssistantInput}
        isBusy={assistantBusy}
        send={sendAssistantChat}
        onNewChat={resetAssistantChat}
      />
      {fileBrowser?.open ? <FileBrowserModal tx={tx} state={fileBrowser} setState={setFileBrowser} api={(path) => api<FileBrowserResponse>(path)} /> : null}
    </main>
  );
}

function SharedConfigPanel(props: {
  tx: Record<string, string>;
  lang: Lang;
  changeLanguage: (value: Lang) => void;
  outputName: string;
  setOutputName: (value: string) => void;
  email: string;
  setEmail: (value: string) => void;
  rememberEmail: boolean;
  setRememberEmail: (value: boolean) => void;
  aiBaseUrl: string;
  setAiBaseUrl: (value: string) => void;
  aiModel: string;
  setAiModel: (value: string) => void;
  offlineModels: OfflineModelsResponse | null;
  detectOfflineModels: () => void;
  aiAutoRun: boolean;
  setAiAutoRun: (value: boolean) => void;
  assayType: "qPCR" | "Conventional";
  setAssayType: (value: "qPCR" | "Conventional") => void;
  parameters: Parameters;
  setParam: (key: string, value: string | number | boolean) => void;
  restore: () => void;
}) {
  return (
    <section className="config-panel">
      <div className="panel-heading">
        <h3>{props.tx.sharedConfig}</h3>
        <button className="ghost" type="button" onClick={props.restore}>
          {props.tx.restore}
        </button>
      </div>

      <div className="config-widget" id="config-general-config">
        <div className="config-widget-body">
          <label>
            {props.tx.language}
            <select value={props.lang} onChange={(event) => props.changeLanguage(event.target.value as Lang)}>
              <option value="vi">Tiếng Việt</option>
              <option value="en">English</option>
            </select>
          </label>
          <label>
            {props.tx.outputName}
            <input value={props.outputName} onChange={(event) => props.setOutputName(event.target.value)} />
          </label>
          <label>
            {props.tx.email}
            <input value={props.email} onChange={(event) => props.setEmail(event.target.value)} />
          </label>
          <label className="toggle-row">
            <span>{props.tx.rememberEmail}</span>
            <input type="checkbox" checked={props.rememberEmail} onChange={(event) => props.setRememberEmail(event.currentTarget.checked)} />
          </label>
          <label>
            {props.tx.assayType}
            <select value={props.assayType} onChange={(event) => props.setAssayType(event.target.value as "qPCR" | "Conventional")}>
              <option value="qPCR">qPCR</option>
              <option value="Conventional">{props.tx.conventional}</option>
            </select>
          </label>
        </div>
      </div>

      <div className="config-widget" id="config-ai-config">
        <div className="config-widget-head">{props.tx.aiConfig}</div>
        <div className="config-widget-body">
          <label>
            {props.tx.aiEndpoint}
            <input value={props.aiBaseUrl} onChange={(event) => props.setAiBaseUrl(event.target.value)} />
          </label>
          <label>
            {props.tx.aiModel}
            {props.offlineModels?.models.length ? (
              <select value={props.aiModel} onChange={(event) => props.setAiModel(event.target.value)}>
                {props.offlineModels.models.map((model) => (
                  <option key={model} value={model}>
                    {model}
                  </option>
                ))}
              </select>
            ) : (
              <input value={props.aiModel} onChange={(event) => props.setAiModel(event.target.value)} />
            )}
          </label>
          <div className="model-status">
            <span>
              {props.tx.modelProvider}: {props.offlineModels?.available ? props.offlineModels.provider : props.tx.modelServerOffline}
            </span>
            <button className="ghost" type="button" onClick={() => void props.detectOfflineModels()}>
              {props.tx.detectModels}
            </button>
          </div>
          <label className="toggle-row">
            <span>{props.tx.aiAutoRun}</span>
            <input type="checkbox" checked={props.aiAutoRun} onChange={(event) => props.setAiAutoRun(event.currentTarget.checked)} />
          </label>
        </div>
      </div>

      {(PARAM_GROUPS as readonly string[]).map((group) => {
        const fields = PARAM_FIELDS.filter((f) => f.group === group);
        if (fields.length === 0) return null;
        return (
          <div className="config-widget" key={group} id={`config-${group}`}>
            <div className="config-widget-head">{props.tx[group] ?? group}</div>
            <div className="config-widget-body">
              {group === "general" && props.parameters.degenerate_primers ? <div className="soft-warning">{props.tx.iupacWarning}</div> : null}
              <div className="config-param-grid">
                {fields.map((field) =>
                  field.type === "boolean" ? (
                    <label className="toggle-row" key={field.key}>
                      <span>{props.tx[field.key] ?? field.key}</span>
                      <input type="checkbox" checked={Boolean(props.parameters[field.key])} onChange={(event) => props.setParam(field.key, event.currentTarget.checked)} />
                    </label>
                  ) : (
                    <label key={field.key}>
                      {props.tx[field.key] ?? field.key}
                      <input
                        type="number"
                        min={field.min}
                        max={field.max}
                        step={field.step}
                        value={String(props.parameters[field.key] ?? "")}
                        onChange={(event) => props.setParam(field.key, Number(event.currentTarget.value))}
                      />
                    </label>
                  ),
                )}
              </div>
            </div>
          </div>
        );
      })}
    </section>
  );
}

function RecentJobs({ tx, jobs, select }: { tx: Record<string, string>; jobs: Job[]; select: (job: Job) => void }) {
  return (
    <section className="recent">
      <div className="section-title">{tx.history}</div>
      {jobs.length === 0 ? (
        <p className="muted">{tx.noJobs}</p>
      ) : (
        jobs.slice(0, 8).map((item) => (
          <button key={item.id} className="job-pill" type="button" onClick={() => select(item)}>
            <span className={`dot ${item.status}`} />
            <span>{tx[item.status] ?? item.status}</span>
            <small>{item.source ? sourceText(tx, item.source) : item.id.slice(-8)}</small>
          </button>
        ))
      )}
    </section>
  );
}

function SystemStatusPanel({ tx, status, refresh }: { tx: Record<string, string>; status: StatusPayload | null; refresh: () => void }) {
  return (
    <section className="recent status-panel">
      <div className="panel-heading compact-heading">
        <div className="section-title">{tx.systemStatus}</div>
        <button className="ghost" type="button" onClick={refresh}>
          {tx.refresh}
        </button>
      </div>
      <div className="status-list">
        <StatusLine label={tx.backend} value={status?.status ?? "-"} />
        <StatusLine label={tx.aiServer} value={status?.ai.available ? status.ai.provider : tx.modelServerOffline} />
        <StatusLine label={tx.jobs} value={String(status?.job_count ?? 0)} />
        <StatusLine label={tx.runStorage} value={formatBytes(status?.run_root_size_bytes ?? 0)} />
        <StatusLine label={tx.diskFree} value={formatBytes(status?.disk.free ?? 0)} />
        <StatusLine label={tx.runRoot} value={status?.run_root ?? "-"} />
        <StatusLine label={tx.fileRoots} value={(status?.file_browser_roots ?? []).join(", ") || "-"} />
      </div>
    </section>
  );
}

function StatusLine({ label, value }: { label: string; value: string }) {
  return (
    <div className="status-line">
      <span>{label}:&nbsp;</span>
      <strong title={value}>{value}</strong>
    </div>
  );
}

function WelcomeCard({ tx }: { tx: Record<string, string> }) {
  return (
    <div className="welcome-card">
      <h3>{tx.aboutTitle}</h3>
      <p>{tx.aboutIntro}</p>
      <div className="workflow-steps">
        <div><strong>1.</strong> {tx.aboutWorkflowDesign}</div>
        <div><strong>2.</strong> {tx.aboutWorkflowValidate}</div>
        <div><strong>3.</strong> {tx.aboutWorkflowReview}</div>
      </div>
      <p className="muted">
        {tx.language === "vi"
          ? "Chọn một chức năng ở thanh điều hướng phía trên để bắt đầu."
          : "Select a workflow from the navigation bar above to get started."}
      </p>
    </div>
  );
}

function DashboardHome(props: { tx: Record<string, string>; jobs: Job[]; selectedJob: Job | null; select: (job: Job) => void }) {
  return (
    <section className="panel">
      <PanelTitle title={props.tx.dashboardTitle} />
      {props.jobs.length === 0 ? (
        <>
          <WelcomeCard tx={props.tx} />
          <div className="empty-state">{props.tx.noJobs}</div>
        </>
      ) : (
        <div className="history-list">
          {props.jobs.map((item) => (
            <button key={item.id} type="button" className={props.selectedJob?.id === item.id ? "history-row active" : "history-row"} onClick={() => props.select(item)}>
              <span className={`dot ${item.status}`} />
              <strong>{item.source ? sourceText(props.tx, item.source) : props.tx.jobLabel}</strong>
              <span>{props.tx[item.status] ?? item.status}</span>
              <code>{item.output_dir}</code>
            </button>
          ))}
        </div>
      )}
    </section>
  );
}

function AboutPage({ tx }: { tx: Record<string, string> }) {
  return (
    <section className="about-page">
      <div className="about-hero panel">
        <div>
          <div className="eyebrow">{tx.about}</div>
          <h3>{tx.aboutTitle}</h3>
          <p>{tx.aboutIntro}</p>
        </div>
      </div>
      <div className="about-grid">
        <article className="panel about-card">
          <span>01</span>
          <h3>{tx.aboutDesignTitle}</h3>
          <p>{tx.aboutDesignText}</p>
        </article>
        <article className="panel about-card">
          <span>02</span>
          <h3>{tx.aboutAiTitle}</h3>
          <p>{tx.aboutAiText}</p>
        </article>
        <article className="panel about-card">
          <span>03</span>
          <h3>{tx.aboutLocalTitle}</h3>
          <p>{tx.aboutLocalText}</p>
        </article>
      </div>
      <section className="panel about-workflow">
        <PanelTitle title={tx.aboutWorkflowTitle} />
        <div className="workflow-steps">
          <div>{tx.aboutWorkflowDesign}</div>
          <div>{tx.aboutWorkflowValidate}</div>
          <div>{tx.aboutWorkflowReview}</div>
        </div>
      </section>
    </section>
  );
}

function LocalMode(props: {
  tx: Record<string, string>;
  isElectron: boolean;
  target: string;
  background: string;
  setTarget: (value: string) => void;
  setBackground: (value: string) => void;
  outputName: string;
  setOutputName: (value: string) => void;
  requestBrowse: (path: string, kind: BrowseKind, onSelect: (path: string) => void) => void;
  parameters: Parameters;
  onSubmit: (event: FormEvent<HTMLFormElement>) => void;
  onUploadSubmit: (event: FormEvent<HTMLFormElement>) => void;
  isBusy: boolean;
  uploadTargetFiles: File[];
  uploadBackgroundFiles: File[];
  setUploadTargetFiles: (files: File[]) => void;
  setUploadBackgroundFiles: (files: File[]) => void;
}) {
  const [dragTarget, setDragTarget] = useState(false);
  const [dragBackground, setDragBackground] = useState(false);
  const targetFileInput = useRef<HTMLInputElement>(null);
  const backgroundFileInput = useRef<HTMLInputElement>(null);

  function handleTargetFiles(files: FileList | null) {
    if (!files) return;
    const fastaFiles = Array.from(files).filter((f) => f.name.endsWith(".fasta") || f.name.endsWith(".fa") || f.name.endsWith(".fna"));
    if (fastaFiles.length === 0) return;
    props.setUploadTargetFiles(fastaFiles);
  }

  function handleBackgroundFiles(files: FileList | null) {
    if (!files) return;
    const fastaFiles = Array.from(files).filter((f) => f.name.endsWith(".fasta") || f.name.endsWith(".fa") || f.name.endsWith(".fna"));
    if (fastaFiles.length === 0) return;
    props.setUploadBackgroundFiles(fastaFiles);
  }

  if (props.isElectron) {
    return (
      <form className="panel" onSubmit={props.onSubmit}>
        <PanelTitle title={props.tx.localTitle} />
        <div className="field-grid">
          <PathInput label={props.tx.targetPath} value={props.target} onChange={props.setTarget} tx={props.tx} requestBrowse={props.requestBrowse} kind="directory" />
          <PathInput label={props.tx.backgroundPath} value={props.background} onChange={props.setBackground} tx={props.tx} requestBrowse={props.requestBrowse} kind="directory" />
        </div>
        <PreflightChecklist
          tx={props.tx}
          items={[
            { label: props.tx.targetPath, ok: Boolean(props.target.trim()) },
            { label: props.tx.backgroundPath, ok: Boolean(props.background.trim()) },
            { label: props.tx.outputName, ok: Boolean(props.outputName.trim()) },
            { label: props.tx.degenerate_primers, ok: true, note: parameterSwitchText(props.parameters.degenerate_primers, props.tx) },
          ]}
        />
        <button className="primary" disabled={props.isBusy} type="submit">
          {props.tx.runLocal}
        </button>
      </form>
    );
  }

  return (
    <form className="panel" onSubmit={props.onUploadSubmit}>
      <PanelTitle title={props.tx.source_upload} />
      <div className="field-grid">
        <div
          className={`drop-zone ${dragTarget ? "drag-over" : ""}`}
          onDragOver={(e) => { e.preventDefault(); setDragTarget(true); }}
          onDragLeave={() => setDragTarget(false)}
          onDrop={(e) => { e.preventDefault(); setDragTarget(false); handleTargetFiles(e.dataTransfer.files); }}
          onClick={() => targetFileInput.current?.click()}
        >
          <input
            ref={targetFileInput}
            type="file"
            multiple
            accept=".fasta,.fa,.fna"
            style={{ display: "none" }}
            onChange={(e) => handleTargetFiles(e.target.files)}
          />
          {props.uploadTargetFiles.length > 0 ? (
            <div>
              <strong>{props.tx.targetPath}:</strong>
              {props.uploadTargetFiles.map((f) => (
                <div key={f.name} className="file-chip">{f.name} ({(f.size / 1024).toFixed(1)} KB)</div>
              ))}
            </div>
          ) : (
            <span className="drop-zone-hint">{props.tx.targetPath} — {props.tx.browse}</span>
          )}
        </div>
        <div
          className={`drop-zone ${dragBackground ? "drag-over" : ""}`}
          onDragOver={(e) => { e.preventDefault(); setDragBackground(true); }}
          onDragLeave={() => setDragBackground(false)}
          onDrop={(e) => { e.preventDefault(); setDragBackground(false); handleBackgroundFiles(e.dataTransfer.files); }}
          onClick={() => backgroundFileInput.current?.click()}
        >
          <input
            ref={backgroundFileInput}
            type="file"
            multiple
            accept=".fasta,.fa,.fna"
            style={{ display: "none" }}
            onChange={(e) => handleBackgroundFiles(e.target.files)}
          />
          {props.uploadBackgroundFiles.length > 0 ? (
            <div>
              <strong>{props.tx.backgroundPath}:</strong>
              {props.uploadBackgroundFiles.map((f) => (
                <div key={f.name} className="file-chip">{f.name} ({(f.size / 1024).toFixed(1)} KB)</div>
              ))}
            </div>
          ) : (
            <span className="drop-zone-hint">{props.tx.backgroundPath} — {props.tx.browse}</span>
          )}
        </div>
      </div>
      <PreflightChecklist
        tx={props.tx}
        items={[
          { label: props.tx.targetPath, ok: props.uploadTargetFiles.length > 0 },
          { label: props.tx.backgroundPath, ok: props.uploadBackgroundFiles.length > 0 },
          { label: props.tx.outputName, ok: Boolean(props.outputName.trim()) },
          { label: props.tx.degenerate_primers, ok: true, note: parameterSwitchText(props.parameters.degenerate_primers, props.tx) },
        ]}
      />
      <button className="primary" disabled={props.isBusy} type="submit">
        {props.tx.runLocal}
      </button>
    </form>
  );
}



function AutoMode(props: {
  tx: Record<string, string>;
  email: string;
  setEmail: (value: string) => void;
  targets: QueryItem[];
  backgrounds: QueryItem[];
  setTargets: (items: QueryItem[]) => void;
  setBackgrounds: (items: QueryItem[]) => void;
  outputName: string;
  setOutputName: (value: string) => void;
  parameters: Parameters;
  onSubmit: (event: FormEvent<HTMLFormElement>) => void;
  isBusy: boolean;
}) {
  return (
    <form className="panel" onSubmit={props.onSubmit}>
      <PanelTitle title={props.tx.autoTitle} />
      <QueryEditor title={props.tx.targetQueries} items={props.targets} setItems={props.setTargets} addLabel={props.tx.addTarget} tx={props.tx} defaultCount={500} />
      <QueryEditor title={props.tx.backgroundQueries} items={props.backgrounds} setItems={props.setBackgrounds} addLabel={props.tx.addBackground} tx={props.tx} defaultCount={50} />
      <PreflightChecklist
        tx={props.tx}
        items={[
          { label: props.tx.email, ok: validEmail(props.email) },
          { label: props.tx.targetQueries, ok: nonEmptyQueries(props.targets).length > 0 },
          { label: props.tx.backgroundQueries, ok: nonEmptyQueries(props.backgrounds).length > 0 },
          { label: props.tx.degenerate_primers, ok: true, note: parameterSwitchText(props.parameters.degenerate_primers, props.tx) },
        ]}
      />
      <button className="primary" disabled={props.isBusy} type="submit">
        {props.tx.runAuto}
      </button>
    </form>
  );
}

function AiMode(props: {
  tx: Record<string, string>;
  lang: Lang;
  email: string;
  setEmail: (value: string) => void;
  aiBaseUrl: string;
  setAiBaseUrl: (value: string) => void;
  aiModel: string;
  setAiModel: (value: string) => void;
  aiAutoRun: boolean;
  setAiAutoRun: (value: boolean) => void;
  chat: ChatMessage[];
  chatInput: string;
  setChatInput: (value: string) => void;
  sendChat: (event: FormEvent<HTMLFormElement>) => void;
  proposal: Record<string, unknown> | null;
  blockedReason: string;
  runProposal: () => void;
  isBusy: boolean;
  onNewChat: () => void;
}) {
  return (
    <section className="panel ai-panel">
      <div className="panel-heading ai-heading">
        <div className="action-group">
          <button className="ghost" type="button" title={props.tx.newChat} onClick={props.onNewChat}>
            ✚
          </button>
          <button className="ghost" type="button" title={props.tx.copyAll} onClick={() => void copyText(chatTranscript(props.chat, props.tx))}>
            ⎘
          </button>
        </div>
      </div>
      <div className="chat-box">
        {props.chat.map((message, index) => (
          <div key={index} className={`chat-message ${message.role}`}>
            <strong>{message.role === "user" ? props.tx.userLabel : props.tx.aiLabel}</strong>
            <p>{renderChatContent(message.content)}</p>
          </div>
        ))}
        {props.isBusy ? <ThinkingIndicator label={props.tx.thinking} /> : null}
      </div>

      <form className="chat-input" onSubmit={props.sendChat}>
        <input value={props.chatInput} placeholder={props.tx.chatPlaceholder} onChange={(event) => props.setChatInput(event.target.value)} />
        <button className="primary" disabled={props.isBusy} type="submit">
          {props.tx.send}
        </button>
      </form>

      {props.blockedReason ? <div className="alert">{props.tx.blocked}: {props.blockedReason}</div> : null}
      {props.proposal ? (
        <div className="proposal">
          <div className="panel-heading">
            <h3>{props.tx.proposal}</h3>
            <button className="primary" type="button" onClick={props.runProposal}>
              {props.tx.runProposal}
            </button>
          </div>
          <ProposalSummary proposal={props.proposal} tx={props.tx} />
        </div>
      ) : null}
    </section>
  );
}

function FloatingAssistant(props: {
  tx: Record<string, string>;
  open: boolean;
  setOpen: (value: boolean) => void;
  messages: ChatMessage[];
  input: string;
  setInput: (value: string) => void;
  isBusy: boolean;
  send: (event: FormEvent<HTMLFormElement>) => void;
  onNewChat: () => void;
}) {
  return (
    <aside className={props.open ? "assistant-dock open" : "assistant-dock"} aria-label={props.tx.assistantDock}>
      <button className="assistant-tab" type="button" onClick={() => props.setOpen(!props.open)} aria-expanded={props.open}>
        <span>AI</span>
        {props.open ? props.tx.assistantClose : props.tx.assistantOpen}
      </button>
      {props.open ? (
        <section className="assistant-card">
          <div className="assistant-header">
            <div>
              <strong>{props.tx.assistantDock}</strong>
              <span>{props.tx.aiConfig}</span>
            </div>
            <div className="action-group">
              <button className="ghost" type="button" title={props.tx.newChat} onClick={props.onNewChat}>
                ✚
              </button>
              <button className="ghost" type="button" title={props.tx.copyAll} onClick={() => void copyText(chatTranscript(props.messages, props.tx))}>
                ⎘
              </button>
            </div>
          </div>
          <div className="assistant-messages">
            {props.messages.map((message, index) => (
              <div key={index} className={`chat-message ${message.role}`}>
                <strong>{message.role === "user" ? props.tx.userLabel : props.tx.aiLabel}</strong>
                <p>{renderChatContent(message.content)}</p>
              </div>
            ))}
            {props.isBusy ? <ThinkingIndicator label={props.tx.thinking} /> : null}
          </div>
          <form className="chat-input" onSubmit={props.send}>
            <input value={props.input} placeholder={props.tx.assistantPlaceholder} onChange={(event) => props.setInput(event.target.value)} />
            <button className="primary" disabled={props.isBusy} type="submit">
              {props.tx.send}
            </button>
          </form>
        </section>
      ) : null}
    </aside>
  );
}

function ThinkingIndicator({ label }: { label: string }) {
  return (
    <div className="chat-message assistant thinking-message" aria-live="polite">
      <strong>AI</strong>
      <p>
        <span className="thinking-dots" aria-hidden="true">
          <span />
          <span />
          <span />
        </span>
        {label}
      </p>
    </div>
  );
}

function renderChatContent(content: string): ReactNode[] {
  const parts = content.split(/(\*\*[^*]+\*\*)/g);
  return parts.map((part, index) => {
    if (part.startsWith("**") && part.endsWith("**")) {
      return (
        <span key={index} className="chat-emphasis">
          {part.slice(2, -2)}
        </span>
      );
    }
    return part;
  });
}

async function copyText(value: string) {
  await navigator.clipboard?.writeText(value);
}

function chatTranscript(messages: ChatMessage[], tx: Record<string, string>) {
  return messages.map((message) => `${message.role === "user" ? tx.userLabel : tx.aiLabel}: ${message.content}`).join("\n\n");
}

function parameterSwitchText(value: string | number | boolean | undefined, tx: Record<string, string>) {
  return value ? tx.enabledItem : tx.disabledItem;
}

function FileBrowserModal(props: {
  tx: Record<string, string>;
  state: FileBrowserState;
  setState: (state: FileBrowserState | null) => void;
  api: (path: string) => Promise<FileBrowserResponse>;
}) {
  const [payload, setPayload] = useState<FileBrowserResponse | null>(null);
  const [pathInput, setPathInput] = useState(props.state.path);
  const [browserError, setBrowserError] = useState("");

  useEffect(() => {
    let active = true;
    setBrowserError("");
    props
      .api(`/api/files/browse?kind=${props.state.kind}&path=${encodeURIComponent(props.state.path)}`)
      .then((next) => {
        if (!active) return;
        setPayload(next);
        setPathInput(next.path);
      })
      .catch((err) => {
        if (active) setBrowserError(messageOf(err, props.tx.unknownError));
      });
    return () => {
      active = false;
    };
  }, [props.state.path, props.state.kind]);

  function go(path: string) {
    props.setState({ ...props.state, path });
  }

  function choose(path: string) {
    props.state.onSelect(path);
    props.setState(null);
  }

  function submitPath(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    go(pathInput);
  }

  return (
    <div className="modal-backdrop" role="dialog" aria-modal="true" aria-label={props.tx.fileBrowser}>
      <section className="file-browser panel">
        <div className="panel-heading">
          <h3>{props.tx.fileBrowser}</h3>
          <button className="ghost" type="button" onClick={() => props.setState(null)}>
            {props.tx.assistantClose}
          </button>
        </div>
        <form className="browse-location" onSubmit={submitPath}>
          <label>
            {props.tx.currentFolder}
            <input value={pathInput} onChange={(event) => setPathInput(event.target.value)} />
          </label>
          <button className="primary" type="submit">
            {props.tx.openFolder}
          </button>
        </form>
        {browserError ? <div className="alert">{browserError}</div> : null}
        <div className="browser-actions">
          {payload?.parent ? (
            <button className="ghost" type="button" onClick={() => go(payload.parent ?? "")}>
              {props.tx.parentFolder}
            </button>
          ) : null}
          {payload ? (
            <button className="primary" type="button" onClick={() => choose(payload.path)}>
              {props.tx.choose}: {payload.path}
            </button>
          ) : null}
        </div>
        <div className="browser-list">
          {(payload?.entries ?? []).map((entry) => (
            <div key={entry.path} className={entry.is_dir ? "browser-row directory" : "browser-row"}>
              <button className="browser-name" type="button" onClick={() => (entry.is_dir ? go(entry.path) : choose(entry.path))}>
                <span>{entry.is_dir ? props.tx.folder : props.tx.file}</span>
                <strong>{entry.name}</strong>
                <code>{entry.path}</code>
              </button>
              <button className="ghost browse-button" type="button" onClick={() => choose(entry.path)}>
                {props.tx.choose}
              </button>
            </div>
          ))}
        </div>
      </section>
    </div>
  );
}

function ValidationMode(props: {
  tx: Record<string, string>;
  mode: DataMode;
  setMode: (value: DataMode) => void;
  action: "validate" | "probe";
  setAction: (value: "validate" | "probe") => void;
  primers: PrimerPair[];
  setPrimers: (items: PrimerPair[]) => void;
  targetPath: string;
  backgroundPath: string;
  setTargetPath: (value: string) => void;
  setBackgroundPath: (value: string) => void;
  requestBrowse: (path: string, kind: BrowseKind, onSelect: (path: string) => void) => void;
  email: string;
  setEmail: (value: string) => void;
  targets: QueryItem[];
  backgrounds: QueryItem[];
  setTargets: (items: QueryItem[]) => void;
  setBackgrounds: (items: QueryItem[]) => void;
  outputName: string;
  setOutputName: (value: string) => void;
  maxMismatch: number;
  setMaxMismatch: (value: number) => void;
  maxLen: number;
  setMaxLen: (value: number) => void;
  extractSeq: boolean;
  setExtractSeq: (value: boolean) => void;
  parameters: Parameters;
  onSubmit: (event: FormEvent<HTMLFormElement>) => void;
  isBusy: boolean;
}) {
  return (
    <form className="panel" onSubmit={props.onSubmit}>
      <PanelTitle title={props.tx.validateTitle} />
      <div style={{ display: "flex", flexDirection: "column", gap: "0.75rem" }}>
        <Segmented<"validate" | "probe">
          value={props.action}
          options={["validate", "probe"]}
          label={props.tx.validationAction}
          tx={{
            ...props.tx,
            validate: props.tx.validateMode,
            probe: props.tx.probeMode,
          }}
          onChange={props.setAction}
        />
        <Segmented<DataMode> value={props.mode} options={["local", "auto"]} label={props.tx.sourceMode} tx={props.tx} onChange={props.setMode} />
      </div>
      <PrimerEditor tx={props.tx} primers={props.primers} setPrimers={props.setPrimers} action={props.action} />
      {props.mode === "local" ? (
        <div className="field-grid">
          <PathInput label={props.tx.targetPath} value={props.targetPath} onChange={props.setTargetPath} tx={props.tx} requestBrowse={props.requestBrowse} kind="directory" />
          <PathInput label={props.tx.backgroundPath} value={props.backgroundPath} onChange={props.setBackgroundPath} tx={props.tx} requestBrowse={props.requestBrowse} kind="directory" />
        </div>
      ) : (
        <>
          <QueryEditor title={props.tx.targetQueries} items={props.targets} setItems={props.setTargets} addLabel={props.tx.addTarget} tx={props.tx} defaultCount={50} />
          <QueryEditor title={props.tx.backgroundQueries} items={props.backgrounds} setItems={props.setBackgrounds} addLabel={props.tx.addBackground} tx={props.tx} defaultCount={10} />
        </>
      )}
      <div className="field-grid compact">
        <label>
          {props.tx.maxMismatch}
          <input type="number" min={0} max={12} value={props.maxMismatch} onChange={(event) => props.setMaxMismatch(Number(event.target.value))} />
        </label>
        <label>
          {props.tx.maxLen}
          <input type="number" min={50} max={10000} value={props.maxLen} onChange={(event) => props.setMaxLen(Number(event.target.value))} />
        </label>
        <label className="toggle-row">
          <span>{props.tx.extractSeq}</span>
          <input type="checkbox" checked={props.extractSeq} onChange={(event) => props.setExtractSeq(event.currentTarget.checked)} />
        </label>
      </div>
      <PreflightChecklist
        tx={props.tx}
        items={[
          { label: props.tx.addPrimer, ok: props.primers.some(validPrimerPair) },
          { label: props.mode === "local" ? props.tx.targetPath : props.tx.targetQueries, ok: props.mode === "local" ? Boolean(props.targetPath.trim()) : nonEmptyQueries(props.targets).length > 0 },
          { label: props.mode === "local" ? props.tx.backgroundPath : props.tx.backgroundQueries, ok: props.mode === "local" ? Boolean(props.backgroundPath.trim()) : true },
          { label: props.tx.maxMismatch, ok: props.maxMismatch >= 0 },
        ]}
      />
      <button className="primary" disabled={props.isBusy} type="submit">
        {props.tx.runValidation}
      </button>
    </form>
  );
}

function MultiplexModePanel(props: {
  tx: Record<string, string>;
  mode: MultiplexMode;
  setMode: (value: MultiplexMode) => void;
  assayType: "qPCR" | "Conventional";
  setAssayType: (value: "qPCR" | "Conventional") => void;
  localTargets: string[];
  setLocalTargets: (items: string[]) => void;
  backgroundPath: string;
  setBackgroundPath: (value: string) => void;
  requestBrowse: (path: string, kind: BrowseKind, onSelect: (path: string) => void) => void;
  email: string;
  setEmail: (value: string) => void;
  targets: QueryItem[];
  backgrounds: QueryItem[];
  setTargets: (items: QueryItem[]) => void;
  setBackgrounds: (items: QueryItem[]) => void;
  folders: string[];
  setFolders: (items: string[]) => void;
  outputName: string;
  setOutputName: (value: string) => void;
  parameters: Parameters;
  onSubmit: (event: FormEvent<HTMLFormElement>) => void;
  isBusy: boolean;
}) {
  return (
    <form className="panel" onSubmit={props.onSubmit}>
      <PanelTitle title={props.tx.multiplexTitle} />
      <Segmented<MultiplexMode> value={props.mode} options={["local", "auto", "analyze"]} label={props.tx.sourceMode} tx={props.tx} onChange={props.setMode} />
      {props.mode === "local" ? (
        <>
          <PathList title={props.tx.localTargets} items={props.localTargets} setItems={props.setLocalTargets} addLabel={props.tx.addFolder} tx={props.tx} requestBrowse={props.requestBrowse} />
          <PathInput label={props.tx.backgroundPath} value={props.backgroundPath} onChange={props.setBackgroundPath} tx={props.tx} requestBrowse={props.requestBrowse} kind="directory" />
        </>
      ) : null}
      {props.mode === "auto" ? (
        <>
          <QueryEditor title={props.tx.targetQueries} items={props.targets} setItems={props.setTargets} addLabel={props.tx.addTarget} tx={props.tx} defaultCount={500} />
          <QueryEditor title={props.tx.backgroundQueries} items={props.backgrounds} setItems={props.setBackgrounds} addLabel={props.tx.addBackground} tx={props.tx} defaultCount={50} />
        </>
      ) : null}
      {props.mode === "analyze" ? <PathList title={props.tx.existingFolders} items={props.folders} setItems={props.setFolders} addLabel={props.tx.addFolder} tx={props.tx} requestBrowse={props.requestBrowse} /> : null}
      <PreflightChecklist
        tx={props.tx}
        items={[
          {
            label: props.mode === "local" ? props.tx.localTargets : props.mode === "auto" ? props.tx.targetQueries : props.tx.existingFolders,
            ok:
              props.mode === "local"
                ? props.localTargets.filter((item) => item.trim()).length >= 2
                : props.mode === "auto"
                  ? nonEmptyQueries(props.targets).length >= 2
                  : props.folders.filter((item) => item.trim()).length >= 2,
          },
          { label: props.mode === "auto" ? props.tx.email : props.tx.backgroundPath, ok: props.mode === "auto" ? validEmail(props.email) : props.mode === "analyze" ? true : Boolean(props.backgroundPath.trim()) },
          { label: props.tx.assayType, ok: true, note: props.assayType },
          { label: props.tx.degenerate_primers, ok: true, note: parameterSwitchText(props.parameters.degenerate_primers, props.tx) },
        ]}
      />
      <button className="primary" disabled={props.isBusy} type="submit">
        {props.tx.runMultiplex}
      </button>
    </form>
  );
}

function QueryEditor(props: {
  title: string;
  items: QueryItem[];
  setItems: (items: QueryItem[]) => void;
  addLabel: string;
  defaultCount: number;
  tx: Record<string, string>;
}) {
  function update(index: number, patch: Partial<QueryItem>) {
    props.setItems(props.items.map((item, i) => (i === index ? { ...item, ...patch } : item)));
  }
  return (
    <div className="query-editor">
      <div className="subheading">{props.title}</div>
      {props.items.map((item, index) => (
        <div className="query-row" key={index}>
          <input aria-label={props.tx.query} value={item.query} onChange={(event) => update(index, { query: event.target.value })} />
          <input aria-label={props.tx.size} type="number" min={0} step={0.1} value={item.size} onChange={(event) => update(index, { size: Number(event.target.value) })} />
          <input aria-label={props.tx.count} type="number" min={0} step={1} value={item.count} onChange={(event) => update(index, { count: Number(event.target.value) })} />
          <select aria-label={props.tx.type} value={item.type} onChange={(event) => update(index, { type: event.target.value as QueryItem["type"] })}>
            <option value="genome">{props.tx.genome}</option>
            <option value="gene">{props.tx.gene}</option>
          </select>
          <button className="ghost" type="button" onClick={() => props.setItems(props.items.filter((_, i) => i !== index))}>
            {props.tx.remove}
          </button>
        </div>
      ))}
      <button className="ghost" type="button" onClick={() => props.setItems([...props.items, { query: "", size: 0, count: props.defaultCount, type: "genome" }])}>
        {props.addLabel}
      </button>
    </div>
  );
}

function PrimerEditor(props: { tx: Record<string, string>; primers: PrimerPair[]; setPrimers: (items: PrimerPair[]) => void; action?: "validate" | "probe" }) {
  function update(index: number, patch: Partial<PrimerPair>) {
    props.setPrimers(props.primers.map((item, i) => (i === index ? { ...item, ...patch } : item)));
  }
  return (
    <div className="query-editor">
      <div className="subheading">{props.tx.validate}</div>
      {props.primers.map((item, index) => (
        <div className="primer-row" key={index}>
          <input aria-label={props.tx.primerName} placeholder={props.tx.primerName} value={item.name} onChange={(event) => update(index, { name: event.target.value })} />
          <input aria-label={props.tx.fwd} placeholder={props.tx.fwd} value={item.fwd} onChange={(event) => update(index, { fwd: event.target.value })} />
          <input aria-label={props.tx.rev} placeholder={props.tx.rev} value={item.rev} onChange={(event) => update(index, { rev: event.target.value })} />
          {props.action !== "probe" && (
            <input aria-label={props.tx.probe} placeholder={props.tx.probe} value={item.probe || ""} onChange={(event) => update(index, { probe: event.target.value })} />
          )}
          <button className="ghost" type="button" onClick={() => props.setPrimers(props.primers.filter((_, i) => i !== index))}>
            {props.tx.remove}
          </button>
        </div>
      ))}
      <button className="ghost" type="button" onClick={() => props.setPrimers([...props.primers, { name: `M${props.primers.length + 1}`, fwd: "", rev: "", probe: "" }])}>
        {props.tx.addPrimer}
      </button>
    </div>
  );
}

function PathInput(props: {
  label: string;
  value: string;
  onChange: (value: string) => void;
  tx: Record<string, string>;
  requestBrowse: (path: string, kind: BrowseKind, onSelect: (path: string) => void) => void;
  kind?: BrowseKind;
}) {
  return (
    <label>
      {props.label}
      <div className="browse-row">
        <input value={props.value} onChange={(event) => props.onChange(event.target.value)} />
        <button className="ghost browse-button" type="button" onClick={() => props.requestBrowse(props.value, props.kind ?? "directory", props.onChange)}>
          {props.tx.browse}
        </button>
      </div>
    </label>
  );
}

function PathList(props: {
  title: string;
  items: string[];
  setItems: (items: string[]) => void;
  addLabel: string;
  tx: Record<string, string>;
  requestBrowse: (path: string, kind: BrowseKind, onSelect: (path: string) => void) => void;
}) {
  function update(index: number, value: string) {
    props.setItems(props.items.map((item, i) => (i === index ? value : item)));
  }
  return (
    <div className="query-editor">
      <div className="subheading">{props.title}</div>
      {props.items.map((item, index) => (
        <div className="path-row" key={index}>
          <input value={item} onChange={(event) => update(index, event.target.value)} />
          <button className="ghost browse-button" type="button" onClick={() => props.requestBrowse(item, "directory", (path) => update(index, path))}>
            {props.tx.browse}
          </button>
          <button className="ghost" type="button" onClick={() => props.setItems(props.items.filter((_, i) => i !== index))}>
            {props.tx.remove}
          </button>
        </div>
      ))}
      <button className="ghost" type="button" onClick={() => props.setItems([...props.items, ""])}>
        {props.addLabel}
      </button>
    </div>
  );
}

function Segmented<T extends string>(props: { value: T; options: T[]; label: string; tx: Record<string, string>; onChange: (value: T) => void }) {
  return (
    <div className="segmented-row">
      <span>{props.label}</span>
      <div className="segmented">
        {props.options.map((option) => (
          <button key={option} type="button" className={props.value === option ? "selected" : ""} onClick={() => props.onChange(option)}>
            {props.tx[option] ?? option}
          </button>
        ))}
      </div>
    </div>
  );
}

function ParametersPanel(props: {
  tx: Record<string, string>;
  parameters: Parameters;
  setParam: (key: string, value: string | number | boolean) => void;
  restore: () => void;
}) {
  return (
    <section className="panel">
      <div className="panel-heading">
        <h3>{props.tx.params}</h3>
        <button className="ghost" type="button" onClick={props.restore}>
          {props.tx.restore}
        </button>
      </div>
      <div className="param-grid">
        {PARAM_FIELDS.map((field) =>
          field.type === "boolean" ? (
            <label className="toggle-row" key={field.key}>
              <span>{props.tx[field.key] ?? field.key}</span>
              <input type="checkbox" checked={Boolean(props.parameters[field.key])} onChange={(event) => props.setParam(field.key, event.currentTarget.checked)} />
            </label>
          ) : (
            <label key={field.key}>
              {props.tx[field.key] ?? field.key}
              <input
                type="number"
                min={field.min}
                max={field.max}
                step={field.step}
                value={String(props.parameters[field.key] ?? "")}
                onChange={(event) => props.setParam(field.key, Number(event.currentTarget.value))}
              />
            </label>
          ),
        )}
      </div>
    </section>
  );
}

function pipelineSteps(terminal: string): Array<{ label: string; done: boolean }> {
  const lines = terminal.split("\n");
  const stepOrder = [
    { label: "Building datasets", pattern: /STAGE 1|BUILDING DATASETS/ },
    { label: "Designing primers", pattern: /STAGE 2|DESIGNING PRIMERS/ },
    { label: "Filtering candidates", pattern: /STAGE 3|FILTERING/ },
    { label: "Saving assay results", pattern: /STAGE 4|SAVING.*ASSAY/ },
    { label: "Validating in-silico PCR", pattern: /STAGE 5|PCR|VALIDAT/ },
    { label: "Saving validation report", pattern: /STAGE 6|SAVING VALIDATION/ },
  ];
  const seen = new Set<string>();
  const steps: Array<{ label: string; done: boolean }> = [];
  let inPipeline = false;
  for (const line of lines) {
    for (const step of stepOrder) {
      if (step.pattern.test(line) && !seen.has(step.label)) {
        seen.add(step.label);
        if (step.label === "Building datasets") inPipeline = true;
        steps.push({ label: step.label, done: false });
      }
    }
    if (/✅/.test(line) && steps.length > 0) steps[steps.length - 1].done = true;
    if (/❌ FAILED/.test(line) && steps.length > 0) steps[steps.length - 1].done = false;
  }
  if (inPipeline && steps.length > 0) steps[0].done = true;
  return steps;
}

function pipelineError(terminal: string): string | null {
  const lines = terminal.split("\n");
  for (const line of lines) {
    if (/❌|FAILED|⚠️|Could not find/.test(line)) return line.trim();
    if (/Traceback/.test(line)) {
      const idx = lines.indexOf(line);
      return lines.slice(idx, Math.min(idx + 4, lines.length)).join(" ").trim();
    }
  }
  return null;
}

function MonitorPanel(props: { tx: Record<string, string>; job: Job | null; terminal: string; cancel: () => void; deleteJob: () => void; refresh: () => void }) {
  const steps = useMemo(() => pipelineSteps(props.terminal), [props.terminal]);
  const errorMsg = useMemo(() => pipelineError(props.terminal), [props.terminal]);
  const failed = props.job?.status === "failed";
  const terminalRef = useRef<HTMLPreElement>(null);
  function scrollToBottom() {
    terminalRef.current?.scrollTo({ top: terminalRef.current.scrollHeight, behavior: "smooth" });
  }
  useEffect(() => {
    if (props.job?.status === "running") scrollToBottom();
  }, [props.terminal]);
  return (
    <section className="panel">
      <div className="panel-heading">
        <h3>{props.tx.monitor}</h3>
        <div className="actions">
          {props.job?.status === "running" ? (
            <button className="danger" type="button" onClick={props.cancel}>
              {props.tx.stop}
            </button>
          ) : null}
          <button className="ghost" type="button" onClick={() => void copyText(props.terminal)}>
            {props.tx.copyAll}
          </button>
          <button className="ghost" type="button" onClick={props.refresh}>
            {props.tx.refresh}
          </button>
          {props.job ? (
            <a className="ghost button-link" href={`${API_BASE}/api/jobs/${props.job.id}/archive`} target="_blank" rel="noreferrer">
              {props.tx.archive}
            </a>
          ) : null}
          {props.job ? (
            <button className="danger" type="button" onClick={props.deleteJob}>
              {props.tx.deleteJob}
            </button>
          ) : null}
        </div>
      </div>
      {failed && errorMsg ? (
        <div className="error-banner">
          <strong>{props.tx.failed}:</strong> {errorMsg}
        </div>
      ) : null}
      <div className="run-summary">
        <Metric label={props.tx.jobLabel} value={props.job ? props.job.id.slice(-12) : "-"} />
        <Metric label={props.tx.sourceLabel} value={props.job?.source ? sourceText(props.tx, props.job.source) : "-"} />
        <Metric label={props.tx.outputLabel} value={props.job?.output_dir ?? "-"} />
      </div>
      {steps.length > 0 ? (
        <div className="pipeline-steps">
          {steps.map((step, idx) => (
            <div key={`${step.label}-${idx}`} className={`pipeline-step ${step.done ? "done" : props.job?.status === "running" ? "active" : ""}`}>
              <span className="step-indicator">{step.done ? "✓" : props.job?.status === "running" ? "○" : "○"}</span>
              <span>{step.label}</span>
            </div>
          ))}
        </div>
      ) : null}
      <div className="terminal-wrap">
        <pre className="terminal" ref={terminalRef}>{props.terminal || props.tx.logWaiting}</pre>
        {props.terminal && props.terminal.length > 200 ? (
          <button className="terminal-scroll-btn ghost" type="button" onClick={scrollToBottom} title="Scroll to bottom">
            ↓
          </button>
        ) : null}
      </div>
    </section>
  );
}

function ResultsPanel({
  tx,
  results,
  requestBrowse,
}: {
  tx: Record<string, string>;
  results: ResultPayload | null;
  requestBrowse: (path: string, kind: BrowseKind, onSelect: (path: string) => void) => void;
}) {
  return (
    <section className="panel">
      <PanelTitle title={tx.results} />
      <ResultInterpretationCard tx={tx} results={results} />
      <ResultTable preview={results?.final_assays} emptyText={tx.noFinal} />
      <ResultTable preview={results?.known_primer_validation} emptyText={tx.noValidation} />
      <ResultTable preview={results?.multiplex_kits} emptyText={tx.noMultiplex} />
      <div className="files">
        {(results?.files ?? []).map((file) => (
          <div key={file.path} className="file-row">
            <span>{file.name}</span>
            <code>{file.path}</code>
            <button className="ghost browse-button" type="button" onClick={() => requestBrowse(file.path, "file", () => undefined)}>
              {tx.browse}
            </button>
          </div>
        ))}
      </div>
    </section>
  );
}

function PreflightChecklist({ tx, items }: { tx: Record<string, string>; items: Array<{ label: string; ok: boolean; note?: string }> }) {
  return (
    <div className="preflight">
      <div className="subheading">{tx.preflight}</div>
      <div className="preflight-list">
        {items.map((item) => (
          <div key={item.label} className={item.ok ? "preflight-row ok" : "preflight-row missing"}>
            <strong>
              <span>{item.ok ? tx.readyItem : tx.missingItem}: </span>
              {item.label}
            </strong>
            {item.note ? <small> - {item.note}</small> : null}
          </div>
        ))}
      </div>
    </div>
  );
}

function ResultInterpretationCard({ tx, results }: { tx: Record<string, string>; results: ResultPayload | null }) {
  const row = resultPreviewRow(results);
  if (!row) return null;
  const candidate = firstAvailable(row, ["Assay", "Assay_ID", "Marker", "Name", "Kit", "Primer", "Target"]) ?? "-";
  const sensitivity = firstAvailable(row, ["Sensitivity", "sensitivity", "Coverage", "Target_Coverage"]) ?? "-";
  const specificity = firstAvailable(row, ["Specificity", "specificity", "Cross_Reactivity", "Non_Target_Hits"]) ?? "-";
  const amplicon = firstAvailable(row, ["Amplicon", "Amplicon_Size", "Product_Size", "product_size"]) ?? "-";
  const primerText = primerCopyText(row);
  return (
    <div className="interpretation-card">
      <div className="panel-heading">
        <h3>{tx.resultInterpretation}</h3>
        {primerText ? (
          <button className="ghost" type="button" onClick={() => void navigator.clipboard?.writeText(primerText)}>
            {tx.copyPrimers}
          </button>
        ) : null}
      </div>
      <div className="interpretation-grid">
        <Metric label={tx.bestCandidate} value={candidate} />
        <Metric label={tx.sensitivity} value={sensitivity} />
        <Metric label={tx.specificity} value={specificity} />
        <Metric label={tx.amplicon} value={amplicon} />
      </div>
    </div>
  );
}

function PanelTitle({ title }: { title: string }) {
  return (
    <div className="panel-title">
      <h3>{title}</h3>
    </div>
  );
}

function Metric({ label, value, icon }: { label: string; value: string; icon?: string }) {
  return (
    <div className="metric">
      <span>{icon ? <span className="metric-icon" aria-hidden="true">{icon}</span> : null}{label}</span>
      <strong title={value}>{value || "-"}</strong>
    </div>
  );
}

function ResultTable({ preview, emptyText }: { preview?: CsvPreview; emptyText: string }) {
  if (!preview?.exists) return <div className="empty-state">{emptyText}</div>;
  if (preview.error) return <div className="empty-state">{preview.error}</div>;
  const columns = preview.columns.slice(0, 8);
  return (
    <div className="table-wrap">
      <table>
        <thead>
          <tr>{columns.map((column) => <th key={column}>{column}</th>)}</tr>
        </thead>
        <tbody>
          {preview.rows.map((row, index) => (
            <tr key={index}>
              {columns.map((column) => <td key={column}>{String(row[column] ?? "")}</td>)}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function quickAiPrompts(tx: Record<string, string>, lang: Lang) {
  if (lang === "en") {
    return [
      {
        label: tx.askLatestResult,
        prompt: "Review the latest design result: select the best assay, identify cross-reactivity risks, and recommend wet-lab next steps.",
      },
      {
        label: tx.askWetlab,
        prompt: "Propose a wet-lab optimization plan for this qPCR assay: primer/probe concentrations, cycling program, negative/positive controls, and acceptance criteria.",
      },
      {
        label: tx.askDegenerate,
        prompt: "Design IUPAC degenerate primers for a SNP-rich or multi-strain target, prioritizing coverage while screening backgrounds strictly.",
      },
      {
        label: tx.askMultiplex,
        prompt: "Create a multiplex PCR/qPCR strategy: choose targets, backgrounds, cross-dimer checks, and balanced Tm constraints.",
      },
    ];
  }
  return [
    {
      label: tx.askLatestResult,
      prompt: "Hãy đánh giá kết quả thiết kế mới nhất: chọn assay tốt nhất, rủi ro cross-reactivity, và khuyến nghị wet-lab.",
    },
    {
      label: tx.askWetlab,
      prompt: "Hãy đề xuất kế hoạch tối ưu wet-lab cho assay qPCR này: nồng độ primer/probe, chương trình nhiệt, kiểm soát âm/dương và tiêu chí nghiệm thu.",
    },
    {
      label: tx.askDegenerate,
      prompt: "Hãy thiết kế mồi thoái hóa IUPAC cho target có biến dị SNP/đa chủng, ưu tiên độ bao phủ nhưng vẫn kiểm tra background chặt chẽ.",
    },
    {
      label: tx.askMultiplex,
      prompt: "Hãy tạo chiến lược multiplex PCR/qPCR: chọn target, background, kiểm tra dimer chéo và cân bằng Tm.",
    },
  ];
}

function resultPreviewRow(results: ResultPayload | null) {
  return results?.final_assays?.rows[0] ?? results?.known_primer_validation?.rows[0] ?? results?.multiplex_kits?.rows[0] ?? null;
}

function firstAvailable(row: Record<string, unknown>, keys: string[]) {
  for (const key of keys) {
    const value = row[key];
    if (value !== undefined && value !== null && String(value).trim()) return String(value);
  }
  return "";
}

function primerCopyText(row: Record<string, unknown>) {
  const fields: Array<[string, string[]]> = [
    ["Forward", ["Forward Primer", "Forward_Primer", "Forward", "Fwd", "fwd"]],
    ["Reverse", ["Reverse Primer", "Reverse_Primer", "Reverse", "Rev", "rev"]],
    ["Probe", ["Probe", "Probe_Seq", "Probe Sequence", "probe"]],
  ];
  return fields
    .map(([label, keys]) => {
      const value = firstAvailable(row, keys);
      return value ? `${label}: ${value}` : "";
    })
    .filter(Boolean)
    .join("\n");
}

function viewTitle(view: WorkspaceView, tx: Record<string, string>) {
  if (view === "dashboard") return tx.dashboardTitle;
  if (view === "local") return tx.localTitle;
  if (view === "auto") return tx.autoTitle;
  if (view === "ai") return tx.aiTitle;
  if (view === "validate") return tx.validateTitle;
  if (view === "about") return tx.aboutTitle;
  return tx.multiplexTitle;
}

function nonEmptyQueries(items: QueryItem[]) {
  return items.filter((item) => item.query.trim());
}

function validEmail(value: string) {
  return /.+@.+\..+/.test(value.trim());
}

function validPrimerPair(item: PrimerPair) {
  const allowed = /^[ATGCRYSWKMNBDHV]+$/i;
  const fwdOk = allowed.test(item.fwd.trim()) && item.fwd.trim().length >= 8;
  const revOk = allowed.test(item.rev.trim()) && item.rev.trim().length >= 8;
  const probeOk = !item.probe || !item.probe.trim() || (allowed.test(item.probe.trim()) && item.probe.trim().length >= 8);
  return fwdOk && revOk && probeOk;
}

function ProposalSummary({ proposal, tx }: { proposal: Record<string, unknown>; tx: Record<string, string> }) {
  const action = String(proposal.action ?? "");
  const items: Array<{ label: string; value: string }> = [];

  if (action === "propose_local_design") {
    items.push({ label: tx.targetPath, value: String(proposal.local_target ?? "") });
    items.push({ label: tx.backgroundPath, value: String(proposal.local_background ?? "") });
  } else if (action === "propose_design" || action === "propose_validation") {
    const targets = proposalList(proposal, "targets", "target");
    const backgrounds = proposalList(proposal, "backgrounds", "background");
    if (targets.length > 0) items.push({ label: tx.targetQueries, value: targets.map((t: unknown) => typeof t === "object" && t ? (t as { query?: string }).query ?? "" : String(t)).join("\n") });
    if (backgrounds.length > 0) items.push({ label: tx.backgroundQueries, value: backgrounds.map((b: unknown) => typeof b === "object" && b ? (b as { query?: string }).query ?? "" : String(b)).join("\n") });
  }

  if (proposal.min_sensitivity) items.push({ label: tx.min_sensitivity, value: `${proposal.min_sensitivity}%` });
  if (proposal.design_min_conservation) items.push({ label: tx.design_min_conservation, value: String(proposal.design_min_conservation) });
  if (proposal.degenerate_primers) items.push({ label: tx.degenerate_primers, value: tx.enabledItem });

  if (items.length === 0) {
    return <pre className="proposal-json">{JSON.stringify(proposal, null, 2)}</pre>;
  }

  return (
    <div className="proposal-items">
      {items.map((item) => (
        <div key={item.label} className="proposal-item">
          <span className="proposal-item-label">{item.label}</span>
          <span className="proposal-item-value">{item.value}</span>
        </div>
      ))}
    </div>
  );
}

function formatBytes(value: number) {
  if (!Number.isFinite(value) || value <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  const index = Math.min(Math.floor(Math.log(value) / Math.log(1024)), units.length - 1);
  return `${(value / 1024 ** index).toFixed(index === 0 ? 0 : 1)} ${units[index]}`;
}

function proposalParameters(proposal: Record<string, unknown>, parameters: Parameters) {
  const next = { ...parameters };
  for (const key of ["design_target_sampling_size", "design_background_sampling_size", "degenerate_primers"]) {
    if (key in proposal) next[key] = proposal[key] as string | number | boolean;
  }
  return next;
}

function proposalList(proposal: Record<string, unknown>, ...keys: string[]) {
  for (const key of keys) {
    const value = proposal[key];
    if (Array.isArray(value)) return value;
  }
  return [];
}

function sourceText(tx: Record<string, string>, source: string) {
  const key = `source_${source.replace(/-/g, "_")}`;
  return tx[key] ?? source;
}

function messageOf(error: unknown, fallback: string) {
  return error instanceof Error ? error.message : fallback;
}
