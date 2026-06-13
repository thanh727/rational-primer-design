const { app, BrowserWindow, ipcMain, Menu, dialog } = require("electron");
const path = require("path");
const fs = require("fs");
const { startBackend, stopBackend, BACKEND_PORT } = require("./python-backend");
const { startFrontendServer, stopFrontendServer, FRONTEND_PORT } = require("./static-server");

// Handle EPIPE error when stdout/stderr streams are closed unexpectedly (e.g., when concurrently terminates)
process.on("uncaughtException", (err) => {
  if (err.code === "EPIPE") {
    process.exit(0);
  }
  console.error("Uncaught Exception:", err);
});


const isDev = process.env.NODE_ENV === "development" || process.argv.includes("--dev");
const BACKEND_URL = `http://127.0.0.1:${BACKEND_PORT}`;

let mainWindow = null;
let loadingWindow = null;

// Single Instance Lock
const gotTheLock = app.requestSingleInstanceLock();
if (!gotTheLock) {
  app.quit();
} else {
  app.on("second-instance", async () => {
    try {
      console.log("[electron] Second instance detected. Ensuring backend is running...");
      await startBackend();
      
      if (!mainWindow) {
        createMainWindow();
      } else {
        if (mainWindow.isMinimized()) mainWindow.restore();
        mainWindow.focus();
        mainWindow.webContents.send("backend-ready", BACKEND_URL);
        mainWindow.reload();
      }
    } catch (err) {
      console.error("[electron] Failed to restart backend on second instance:", err.message);
    }
  });
}


function createLoadingWindow() {
  loadingWindow = new BrowserWindow({
    width: 400,
    height: 300,
    frame: false,
    transparent: false,
    resizable: false,
    center: true,
    title: "Starting...",
    backgroundColor: "#0a0a0a",
  });

  const html = `<!DOCTYPE html>
<html><body style="display:flex;align-items:center;justify-content:center;height:100vh;margin:0;background:#0f172a;color:#e2e8f0;font-family:system-ui,sans-serif;flex-direction:column;gap:16px">
<div style="font-size:48px">🧬</div>
<div style="font-size:18px;font-weight:600">Rational Primer Design</div>
<div style="font-size:13px;color:#94a3b8">Starting backend services...</div>
<div style="width:120px;height:3px;background:#1e293b;border-radius:2px;overflow:hidden">
<div style="width:30%;height:100%;background:#3b82f6;border-radius:2px;animation:pulse 1.2s ease-in-out infinite"></div>
</div>
<style>@keyframes pulse{0%{transform:translateX(-100%)}100%{transform:translateX(400%)}}</style>
</body></html>`;

  loadingWindow.loadURL(`data:text/html;charset=utf-8,${encodeURIComponent(html)}`);
}

function createMainWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 1024,
    minHeight: 700,
    title: "Rational Primer Design",
    show: false,
    webPreferences: {
      preload: path.join(__dirname, "preload.js"),
      contextIsolation: true,
      nodeIntegration: false,
    },
  });

  function loadContent() {
    if (isDev) {
      mainWindow.loadURL("http://localhost:3001");
    } else {
      mainWindow.loadURL(`http://127.0.0.1:${FRONTEND_PORT}`);
    }
  }

  loadContent();

  mainWindow.webContents.on("did-finish-load", () => {
    mainWindow.webContents.send("backend-ready", BACKEND_URL);
  });

  mainWindow.webContents.on("did-fail-load", () => {
    if (isDev) {
      setTimeout(loadContent, 1000);
    }
  });

  mainWindow.once("ready-to-show", () => {
    if (loadingWindow) {
      loadingWindow.close();
      loadingWindow = null;
    }
    mainWindow.show();
    if (isDev) mainWindow.webContents.openDevTools();
  });

  mainWindow.on("closed", () => {
    mainWindow = null;
  });
}

function buildMenu() {
  const template = [
    {
      label: "Rational Primer Design",
      submenu: [
        { role: "about" },
        { type: "separator" },
        { role: "hide" },
        { role: "hideOthers" },
        { role: "unhide" },
        { type: "separator" },
        { role: "quit" },
      ],
    },
    {
      label: "Edit",
      submenu: [
        { role: "undo" },
        { role: "redo" },
        { type: "separator" },
        { role: "cut" },
        { role: "copy" },
        { role: "paste" },
        { role: "selectAll" },
      ],
    },
    {
      label: "View",
      submenu: [
        { role: "reload" },
        { role: "toggleDevTools" },
        { type: "separator" },
        { role: "resetZoom" },
        { role: "zoomIn" },
        { role: "zoomOut" },
        { type: "separator" },
        { role: "togglefullscreen" },
      ],
    },
    {
      label: "Window",
      submenu: [{ role: "minimize" }, { role: "zoom" }, { role: "close" }],
    },
  ];

  if (process.platform === "darwin") {
    Menu.setApplicationMenu(Menu.buildFromTemplate(template));
  } else {
    Menu.setApplicationMenu(Menu.buildFromTemplate(template.slice(1)));
  }
}

ipcMain.handle("get-backend-url", () => BACKEND_URL);
ipcMain.handle("get-app-version", () => app.getVersion());

if (gotTheLock) {
  app.whenReady().then(async () => {
    buildMenu();
    createLoadingWindow();

    try {
      console.log("[electron] Starting Python backend...");
      await startBackend();
      console.log("[electron] Backend ready at", BACKEND_URL);
    } catch (err) {
      console.error("[electron] Backend failed:", err.message);
      dialog.showErrorBox(
        "Backend Connection Error",
        `Failed to start Python backend services:\n${err.message}\n\nThe application will now close.`
      );
      app.quit();
      return;
    }

    if (!isDev) {
      const frontendDir = path.join(__dirname, "..", "frontend", "out");
      if (fs.existsSync(frontendDir)) {
        try {
          await startFrontendServer(frontendDir);
        } catch (err) {
          console.error("[electron] Static server failed:", err.message);
        }
      } else {
        console.error("[electron] Frontend build not found at", frontendDir);
        console.error("[electron] Run 'npm run build:frontend' first");
      }
    }

    createMainWindow();

    app.on("activate", () => {
      if (BrowserWindow.getAllWindows().length === 0) createMainWindow();
    });
  });

  app.on("window-all-closed", () => {
    if (process.platform !== "darwin") app.quit();
  });

  app.on("before-quit", () => {
    stopBackend();
    stopFrontendServer();
  });
}

