# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is Wanjun Gu's personal portfolio website hosted on GitHub Pages at `https://broccolito.github.io/`. It is a static site with no build step — all files are served directly.

## Deployment

Changes are deployed by pushing to the `main` branch. Two sync scripts automate this:

```bash
bash auto_sync.sh      # macOS/Linux
auto_sync.bat          # Windows
```

Both stage all changes (`git add .`), commit with message "automatic commit", pull from origin (allowing unrelated histories), and push to main.

## Architecture

The site is three pages sharing an identical layout and CSS:

- **`index.html`** — About section with bio text
- **`software.html`** — Numbered list of R packages and tools (tKOI, pKOI, KDPS, ggmugs, bolt4jr, gptr, BioRouter, MedCP), each linking out to external documentation
- **`publications.html`** — Publication list

All pages share the same inline CSS design system (CSS custom properties: `--accent`, `--dark`, `--mid`, `--light`, `--border`, `--bg`) and an identical sidebar (sticky profile photo, nav with active state, external profile links). There is no shared CSS file — styles are duplicated inline in each page.

**Assets:** Package logos at `pictures/`, profile photo at `pictures/wanjun.jpg`, CV at `Wanjun Gu CV.pdf`.

## Key Design Decisions

- No JavaScript framework or build toolchain — pure HTML/CSS with no CDN dependencies (all styles are inline)
- Software documentation is hosted externally (GitHub Pages for individual R package repos via `pkgdown`); only logos are stored locally
- When adding a new page, copy the full `<style>` block and sidebar HTML from an existing page to maintain consistency
