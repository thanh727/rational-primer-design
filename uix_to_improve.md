# UX/UI Review — Rational Primer Design

## Visual Design

### Color System
- **Good**: Teal (`#0f8b7d`) + blue (`#2563eb`) accent palette feels trustworthy and scientific
- **Issue**: Purple (`#7c3aed`) used for floating assistant — clashes with the teal/blue scheme. Should be a derivative of the main palette
- **Issue**: No light/dark mode toggle. Dark terminal is good but the rest of the app lacks a dark variant
- **Issue**: Background gradient is subtle but the config sidebar gradient (teal → purple) introduces the purple prematurely

### Typography
- **Good**: Inter font, clean hierarchy, emoji icons in logs
- **Issue**: Sidebar labels at 13px are hard to read on retina displays
- **Issue**: No icon library — everything is pure text, making the nav and buttons hard to scan quickly
- **Issue**: Metric card values (`strong`) don't have a distinctive font weight difference from their labels

### Spacing & Layout
- **Good**: 8px grid system, consistent panel padding (20px), generous gap (16px)
- **Issue**: Sidebar at 330px is too narrow for the parameter controls — inputs get squished
- **Issue**: Top nav with 7 items including "About" alongside tools, which should be separated
- **Issue**: The 4-metric card bar appears on every non-about page but shows the same generic metrics regardless of context (e.g., "Jobs" count on the validate page)

## Component-Level Issues

### Navigation
- **No icons**: All nav links are text-only, hard to scan
- **"About" is mixed with tools** — should be in a secondary position or as a ? icon
- **Active state underline animation** is nice but barely visible on hover

### Config Sidebar
- **Too dense**: 15+ controls crammed into 330px width
- **"Core parameters" section** has 10 parameters in a 2-column grid — each input is tiny
- **No collapsible sections** — user sees everything at once, overwhelming
- **"Restore defaults" button** is tiny and easy to miss
- **Status panel shows "-" placeholders** on first load before API responds — flashes of unstyled content

### Dashboard
- **Empty state**: Shows "No job history" and WelcomeCard (good) but the empty-state div placeholder is a massive 150px dashed box that feels like an error
- **Job history rows**: No timestamps formatted as relative ("2 hours ago"), no pagination, no search
- **Monitor panel**: Terminal is excellent (dark theme, monospace) but no "scroll to bottom" button
- **Results panel**: Tables lack sticky headers, no horizontal scroll affordance, no column resize

### Design Forms (Local/Auto)
- **Form layout**: Single-column form with the monitor/results panel beside it — on the design page, the user fills in the form on the left but the results panel on the right is empty, creating a dead zone
- **Path fields**: Input + "Browse" button works but the input value text is tiny
- **Query rows (auto mode)**: 5-column grid in a narrow space — each input is barely 80px wide
- **Pre-flight checklist**: Useful but looks like a debugging tool, not a UX pattern

### Validation Page
- **Primer rows**: 4-column grid (name, fwd, rev, delete) — the sequence inputs are too narrow for 20+ bp primers
- **No "Paste from Excel" button** — users have to type one by one
- **Validation parameters** (max mismatch, max amplicon) are below the primer list but above the submit — should be a dedicated section

### AI Chat Page
- **Proposal display**: Raw `<pre>` JSON block — meaningless to non-technical users
- **Chat box**: Max height 62vh with no resize handle
- **Quick prompt buttons**: Good idea but no tooltips explaining what they do
- **No streaming response**: User waits for full response before seeing anything

### Multiplex Page
- **Three-mode segmented control** (local/auto/analyze) is not explained — user doesn't know which to pick
- **"Analyze" mode** asks for result folders but provides no way to browse existing results

### Floating Assistant
- **Purple color** clashes with main teal theme
- **Always in bottom-right** — could overlap important content on small screens
- **No unread notification badge** when collapsed

### File Browser Modal
- **No file/folder icons** — just text rows
- **No file size formatting** (shows raw bytes)
- **No search within current directory**
- **No "Upload" button** even though the API supports file upload

## Micro-Interactions

### Missing Feedback
- **No loading skeleton** — placeholders jump from "-" to data
- **No toast notifications** — errors shown in a static `.alert` div at the top of the workspace
- **No success animation** when a job completes
- **Copy buttons** give no visual feedback (no "Copied!" tooltip)

### Animation
- **Good**: panel hover lift, thinking dots pulse
- **Missing**: Page transitions, form submission loading state, smooth list updates

## Accessibility

### Issues
- **Language toggle** doesn't persist across page reloads (saved but needs refresh?)
- **Color contrast**: Muted text (#536663 on white) may fail WCAG AA at small sizes
- **No `aria-live` region** for job status updates — screen readers won't announce changes
- **Form labels** use the generic `<label>` pattern but there's no `htmlFor` linking — relies on nesting
- **Keyboard navigation**: Tab order follows visual order but the dense sidebar has too many tab stops

## Recommendations (Priority-Ordered)

### P1 — High Impact, Low Effort
1. Add **SVG icons** to top navigation items
2. Improve **metric cards** to show page-contextual labels
3. Replace AI proposal **raw JSON** with a structured summary card
4. Add **sticky table headers** and horizontal scroll indicator
5. Add **loading shimmer/skeleton** for initial data fetch
6. Fix **floating assistant color** to match teal theme

### P2 — Medium Impact
7. Add **collapsible sections** to the config sidebar
8. Add **relative timestamps** to job history
9. Add **"scroll to bottom"** button to terminal
10. Improve **file browser** with icons and formatted sizes
11. Add **toast notification** system instead of static alert
12. Add **copy feedback** ("Copied!" tooltip)

### P3 — Lower Impact / Higher Effort
13. Full **dark mode** support
14. **Responsive sidebar** — collapse to icons on narrow screens
15. **Drag-and-drop primer import** from CSV/Excel
16. **Live search/filter** in job history
17. **Resizable panels** (form vs monitor)
18. **Streaming AI responses**
