"use client";

import { useEffect } from "react";

export default function Heartbeat() {
  useEffect(() => {
    const API_BASE = process.env.NEXT_PUBLIC_API_BASE_URL || "http://127.0.0.1:8000";
    const sendHeartbeat = () => {
      fetch(`${API_BASE}/api/heartbeat`, { method: "POST" }).catch(() => {});
    };
    
    // Send initial heartbeat
    sendHeartbeat();
    
    // Heartbeat every 3 seconds
    const interval = setInterval(sendHeartbeat, 3000);
    return () => clearInterval(interval);
  }, []);

  return null;
}
