const http = require("http");
const fs = require("fs");
const path = require("path");

const FRONTEND_PORT = 3001;

const MIME = {
  ".html": "text/html",
  ".js": "text/javascript",
  ".css": "text/css",
  ".json": "application/json",
  ".png": "image/png",
  ".jpg": "image/jpeg",
  ".svg": "image/svg+xml",
  ".ico": "image/x-icon",
  ".woff2": "font/woff2",
};

function tryRead(filePath) {
  return new Promise((resolve, reject) => {
    fs.readFile(filePath, (err, data) => {
      if (err) return reject(err);
      resolve(data);
    });
  });
}

async function serveFile(rootDir, req, res) {
  const urlPath = req.url === "/" ? "/index.html" : req.url;
  const rawPath = path.join(rootDir, urlPath);

  try {
    const data = await tryRead(rawPath);
    const ext = path.extname(rawPath);
    res.writeHead(200, { "Content-Type": MIME[ext] || "application/octet-stream" });
    res.end(data);
    return;
  } catch (err) {
    if (err.code !== "ENOENT" && err.code !== "EISDIR") {
      console.error("[static-server] error reading", rawPath, err.message);
    }
  }

  const ext = path.extname(urlPath);
  if (!ext) {
    for (const candidate of [rawPath + ".html", path.join(rootDir, urlPath, "index.html")]) {
      try {
        const data = await tryRead(candidate);
        res.writeHead(200, { "Content-Type": "text/html" });
        res.end(data);
        return;
      } catch {}
    }
  }

  try {
    const data = await tryRead(path.join(rootDir, "index.html"));
    res.writeHead(200, { "Content-Type": "text/html" });
    res.end(data);
  } catch {
    res.writeHead(404);
    res.end("404");
  }
}

let server = null;

function startFrontendServer(rootDir) {
  return new Promise((resolve, reject) => {
    server = http.createServer((req, res) => serveFile(rootDir, req, res));
    server.listen(FRONTEND_PORT, "127.0.0.1", () => {
      console.log(`[static-server] Serving ${rootDir} on http://127.0.0.1:${FRONTEND_PORT}`);
      resolve(FRONTEND_PORT);
    });
    server.on("error", reject);
  });
}

function stopFrontendServer() {
  if (server) {
    server.close();
    server = null;
  }
}

module.exports = { startFrontendServer, stopFrontendServer, FRONTEND_PORT };
