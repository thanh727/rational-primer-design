"use client";

import { useCallback, useEffect, useRef } from "react";

function getApiBase(): string {
  if (typeof window !== "undefined" && (window as any).electronAPI) {
    return "http://127.0.0.1:8000";
  }
  return process.env.NEXT_PUBLIC_API_BASE_URL || "http://127.0.0.1:8000";
}

export default function Heartbeat() {
  const intervalRef = useRef<ReturnType<typeof setInterval> | null>(null);

  const sendHeartbeat = useCallback(() => {
    fetch(`${getApiBase()}/api/heartbeat`, { method: "POST" }).catch(() => {});
  }, []);

  useEffect(() => {
    sendHeartbeat();
    intervalRef.current = setInterval(sendHeartbeat, 3000);
    return () => {
      if (intervalRef.current) clearInterval(intervalRef.current);
    };
  }, [sendHeartbeat]);

  return null;
}
