# 🧠 Transferable Protocol: Kiến trúc Ứng dụng AI-Centric (Lấy LLM làm Vi xử lý trung tâm)

Tài liệu này mô tả một **Giao thức chuyển đổi (Transferable Protocol)** giúp các kỹ sư và nhà phát triển phần mềm thiết kế, xây dựng các ứng dụng mà trong đó **Mô hình Ngôn ngữ Lớn (LLMs)** không chỉ đóng vai trò là một tính năng phụ (chat/hỏi đáp), mà hoạt động như một **Bộ Vi xử lý trung tâm (Cognitive CPU)** hoặc **Hệ điều hành cốt lõi (Core OS)** của toàn bộ hệ thống.

---

## 1. Triết lý Kiến trúc (The Philosophy)

Trong phần mềm truyền thống, "Vi xử lý trung tâm" của ứng dụng là **Logic code** (`if/else`, vòng lặp, Regex, thuật toán quy tắc cứng). 
Trong kiến trúc **AI-Centric**, "Vi xử lý trung tâm" là **LLM**. 
- Code truyền thống (Python, JS) lùi về làm nhiệm vụ **I/O (Input/Output)**: thu thập dữ liệu đầu vào (UI, Database, File), chuyển đổi định dạng và kết xuất kết quả.
- LLM đảm nhận việc **phân tích, đưa ra quyết định, phân loại, và kết luận**.

### 🌟 4 Nguyên tắc Cốt lõi (Core Principles)
1. **Model Agnosticism (Không phụ thuộc Model):** Ứng dụng không được "trói chặt" vào bất kỳ một Model cụ thể nào (như GPT-4 hay Gemini). Phải có một lớp trừu tượng (Abstraction Layer) cho phép chuyển đổi dễ dàng giữa Cloud API (Google, OpenAI) và Local Model (Ollama, LM Studio) chỉ bằng 1 dòng cấu hình.
2. **Structured Output (Đầu ra có cấu trúc):** AI đóng vai trò xử lý dữ liệu hệ thống, do đó đầu ra bắt buộc phải là định dạng máy có thể đọc được (như **JSON**), không phải văn bản thuần (plain text).
3. **Stateless Processing (Xử lý phi trạng thái):** Mỗi lần gọi AI là một hàm độc lập (Pure Function). Toàn bộ ngữ cảnh (Context) phải được đóng gói gọn gàng trong 1 Prompt duy nhất (Zero-shot hoặc Few-shot) để tối ưu hoá hiệu năng và tránh rò rỉ bộ nhớ.
4. **Resilience (Khả năng chịu lỗi):** API có thể sập, Local model có thể quá tải. Bắt buộc phải có cơ chế Retry (thử lại) và Exponential Backoff (chờ giãn cách) ở tầng gọi AI.

---

## 2. Mô hình 4 Lớp (The 4-Layer Architecture)

Bất kỳ Package nào xây dựng theo giao thức này đều tuân thủ mô hình 4 lớp sau:

### Lớp 1: UI & Data Ingestion (Giao diện & Thu thập dữ liệu)
- **Nhiệm vụ:** Tương tác với người dùng, nhận đầu vào (File PDF, Hình ảnh, Text, Âm thanh).
- **Công nghệ gợi ý:** Streamlit (cho công cụ nội bộ), Next.js / React (cho Web App).
- **Hành vi:** Tiền xử lý dữ liệu (ví dụ: cắt nhỏ ảnh, nén dung lượng, convert file sang Base64) để làm "thức ăn" chuẩn bị đút cho AI.

### Lớp 2: Orchestration & Prompting (Điều phối & Nhúng Ngữ cảnh)
- **Nhiệm vụ:** Đây là phần "Phần mềm" thực sự. Nơi biến Dữ liệu đầu vào + Quy tắc nghiệp vụ (Business Rules) thành Lệnh (Prompt).
- **Thành phần:**
  - `Prompt Builder`: Cấu trúc hoá câu hỏi (vd: "Hãy đọc ảnh này và trích xuất các trường sau...").
  - `Schema Definition`: Định nghĩa sẵn JSON Schema yêu cầu AI trả về.

### Lớp 3: The Abstraction Layer (Lớp Trừu tượng AI - Trái tim hệ thống)
- **Nhiệm vụ:** Đóng gói giao thức giao tiếp với AI. 
- **Thiết kế (Interface):** Chỉ cần một hàm chuẩn hoá: `generate(prompt, image) -> dict (JSON)`
- **Triển khai thực tế:**
  - Môi trường Doanh nghiệp / Có Internet: Gọi qua SDK của Cloud AI (Gemini / OpenAI API).
  - Môi trường Bảo mật cao / Offline: Gọi qua chuẩn OpenAI-compatible tới các Local Server (Ollama / vLLM / LM Studio) đang chạy trên máy chủ nội bộ.

### Lớp 4: The Cognitive Engine (Bộ máy Nhận thức - LLM)
- Xử lý ngôn ngữ, thị giác máy tính và logic suy luận dựa trên dữ liệu từ Lớp 3 chuyển tới, trả về chuỗi JSON chính xác.

---

## 3. Template Triển khai (Implementation Blueprint)

Khi bắt đầu một Package mới, hãy copy khung thiết kế sau:

### A. Định nghĩa Interface (Python Ví dụ)
Tất cả các Backend phải kế thừa từ một Abstract Base Class (ABC).

```python
from abc import ABC, abstractmethod

class LLMBackend(ABC):
    @abstractmethod
    def generate_json(self, prompt: str, system_instruction: str = "", image=None) -> dict:
        """Nhận vào Prompt (và Ảnh), trả về Dictionary (JSON parse sẵn)."""
        pass
        
    @abstractmethod
    def test_connection(self) -> bool:
        pass
```

### B. Xây dựng Prompt "Ép khuôn" JSON (Strict JSON Prompting)
Để AI trở thành Vi xử lý, nó không được nói "Xin chào, đây là kết quả của bạn:". Nó chỉ được khạc ra JSON.

> **Công thức Prompt Chuẩn:**
> 1. [Vai trò] Bạn là một chuyên gia phân tích dữ liệu...
> 2. [Nhiệm vụ] Hãy đọc thông tin đầu vào sau và phân loại nó.
> 3. [Quy tắc] Trả về CHỈ ĐỊNH DẠNG JSON. KHÔNG ĐƯỢC GIẢI THÍCH THÊM. KHÔNG DÙNG MARKDOWN TUKT.
> 4. [Schema] Cấu trúc bắt buộc: `{"status": "PASS|FAIL", "reason": "...", "confidence": 0.9}`

### C. Cơ chế Retry thông minh (The Wrapper)
```python
import time
import json

def robust_ai_call(backend: LLMBackend, prompt: str, retries=3):
    for attempt in range(retries):
        try:
            raw_text = backend.generate_json(prompt)
            # Dọn dẹp markdown nếu AI lỡ sinh ra
            raw_text = raw_text.replace("```json", "").replace("```", "").strip()
            return json.loads(raw_text)
        except Exception as e:
            if attempt == retries - 1:
                return {"error": str(e), "status": "FAILED"}
            time.sleep(2 ** attempt) # Exponential backoff
```

---

## 4. Khả năng Mở rộng (Use Cases áp dụng Protocol)

Protocol này có thể được copy-paste nguyên si cấu trúc để xây dựng vô số hệ thống:

1. **Hệ thống Quản lý Bệnh án Thông minh (Smart EMR):**
   - *Input:* Ảnh chụp sổ y bạ viết tay, hoặc ghi âm lời bác sĩ.
   - *Prompt:* "Trích xuất Tên bệnh nhân, Chẩn đoán, Toa thuốc thành JSON."
   - *AI:* Trả về JSON. Hệ thống tự động push JSON này vào CSDL Bệnh viện.
2. **Kế toán Tự động (Auto-Billing):**
   - *Input:* PDF/Ảnh Hoá đơn, Hợp đồng nhập vào liên tục.
   - *Prompt:* "Tìm Mã số Thuế, Số tiền, Ngày xuất hoá đơn. Đối chiếu có hợp lệ không."
3. **Phân loại Yêu cầu Hỗ trợ (Customer Support Router):**
   - *Input:* Email hoặc Ticket hỗ trợ dài ngoằng của khách hàng.
   - *Prompt:* "Xác định Cảm xúc (Tức giận/Bình thường), Bộ phận cần chuyển đến (Kỹ thuật/Thanh toán), Mức độ khẩn cấp (Cao/Thấp)."

---

## 5. Lời kết
Sự khác biệt giữa một "Chatbot" và một "Phần mềm AI chuyên nghiệp" nằm ở **Sự khắt khe trong việc kiểm soát luồng dữ liệu (Data Flow)**. Bằng cách nhốt LLM vào bên trong một *Abstraction Layer* và buộc nó giao tiếp bằng *JSON*, bạn đã biến một mô hình sinh ngữ ngẫu nhiên thành một con Chip xử lý vô cùng mạnh mẽ, dễ dàng tích hợp, nâng cấp và bảo trì cho bất kỳ Domain nghiệp vụ nào.
