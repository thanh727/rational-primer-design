const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const http = require("http");

const BACKEND_PORT = 8000;
const HEALTH_PATH = "/api/health";
const POLL_INTERVAL = 300;
const POLL_TIMEOUT = 20000;
const MAX_RETRIES = 3;

let backendProcess = null;
let retryCount = 0;

function resolveBinary() {
  const resources = process.resourcesPath || path.join(__dirname, "..");
  const ext = process.platform === "win32" ? ".exe" : "";
  const prodBinary = path.join(resources, "backend", `api_server${ext}`);

  if (process.env.NODE_ENV === "development" || !fs.existsSync(prodBinary)) {
    const venvPython = path.join(__dirname, "..", "venv", "bin", "python");
    if (fs.existsSync(venvPython)) {
      return { command: venvPython, args: ["-m", "rational_design.api_server"] };
    }
    const venvPythonWin = path.join(__dirname, "..", "venv", "Scripts", "python.exe");
    if (fs.existsSync(venvPythonWin)) {
      return { command: venvPythonWin, args: ["-m", "rational_design.api_server"] };
    }
    return { command: "python3", args: ["-m", "rational_design.api_server"] };
  }
  return { command: prodBinary, args: [] };
}


function probe(url) {
  return new Promise((resolve) => {
    const req = http.get(url, (res) => {
      let body = "";
      res.on("data", (c) => (body += c));
      res.on("end", () => resolve(res.statusCode === 200));
    });
    req.on("error", () => resolve(false));
    req.setTimeout(2000, () => { req.destroy(); resolve(false); });
    req.end();
  });
}

function waitForServer(url, timeout) {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    function poll() {
      if (Date.now() - start > timeout) {
        return reject(new Error("Backend did not start within timeout"));
      }
      probe(url).then((ok) => (ok ? resolve() : setTimeout(poll, POLL_INTERVAL)));
    }
    poll();
  });
}

function startBackend() {
  if (backendProcess) return Promise.resolve();

  const { command, args } = resolveBinary();
  const cwd = path.join(__dirname, "..");

  return new Promise((resolve, reject) => {
    backendProcess = spawn(command, args, {
      cwd,
      detached: true,
      stdio: ["ignore", "pipe", "pipe"],
      env: {
        ...process.env,
        RPD_API_HOST: "127.0.0.1",
        RPD_API_PORT: String(BACKEND_PORT),
        RPD_RUN_ROOT: path.join(cwd, "runs", "web_jobs"),
      },
    });

    const log = (prefix) => (data) => {
      const lines = data.toString().trim().split("\n").filter(Boolean);
      lines.forEach((l) => console.log(`[${prefix}] ${l}`));
    };

    backendProcess.stdout.on("data", log("backend"));
    backendProcess.stderr.on("data", log("backend"));

    backendProcess.on("error", (err) => {
      console.error("[backend] spawn error:", err.message);
      backendProcess = null;
      reject(err);
    });

    backendProcess.on("exit", (code, signal) => {
      console.log(`[backend] exited code=${code} signal=${signal}`);
      backendProcess = null;
    });

    waitForServer(`http://127.0.0.1:${BACKEND_PORT}${HEALTH_PATH}`, POLL_TIMEOUT)
      .then(resolve)
      .catch(reject);
  });
}

function stopBackend() {
  if (!backendProcess) return;
  try {
    if (process.platform === "win32") {
      spawn("taskkill", ["/pid", String(backendProcess.pid), "/f", "/t"]);
    } else {
      process.kill(-backendProcess.pid, "SIGTERM");
    }
  } catch {
    try { backendProcess.kill("SIGTERM"); } catch {}
  }
  backendProcess = null;
}

module.exports = { startBackend, stopBackend, BACKEND_PORT };
