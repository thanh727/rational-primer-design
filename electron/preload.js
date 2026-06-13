const { contextBridge, ipcRenderer } = require("electron");

contextBridge.exposeInMainWorld("electronAPI", {
  getBackendUrl: () => ipcRenderer.invoke("get-backend-url"),
  onBackendReady: (callback) => ipcRenderer.on("backend-ready", (_event, url) => callback(url)),
  getAppVersion: () => ipcRenderer.invoke("get-app-version"),
  platform: process.platform,
});
