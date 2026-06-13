import type { Metadata } from "next";
import type { ReactNode } from "react";
import "./globals.css";
import Heartbeat from "../components/Heartbeat";

export const metadata: Metadata = {
  title: "Rational Primer Design",
  description: "Primer and TaqMan assay design workspace",
};

export default function RootLayout({ children }: Readonly<{ children: ReactNode }>) {
  return (
    <html lang="vi">
      <body>
        <Heartbeat />
        {children}
      </body>
    </html>
  );
}

