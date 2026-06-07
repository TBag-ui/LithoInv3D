// svg_plot.h — Minimal SVG chart generation (public domain / CC0)
//
// Zero dependencies.  Writes SVG text that renders in any browser.
// Usage:
//   SVGPlot plot(800, 600);
//   plot.addLine(x, y, n, "blue", "Data");
//   plot.addScatter(x, y, n, "red", 3.0);
//   plot.addHLine(y, "gray", 1.5, "dashed");
//   plot.save("output.svg");
//
// Rows/cols create subplot grids:  SVGPlot plot(800, 600, 2, 2);
#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iomanip>

namespace litho_invert {

struct SVGPlot {
    struct Panel {
        std::string title;
        std::string xlabel;
        std::string ylabel;
        double xmin = 0, xmax = 1, ymin = 0, ymax = 1;
        bool logY = false;

        struct Line {
            std::vector<double> x, y;
            std::string color;
            std::string label;
            double width = 1.5;
            std::string dash;  // "dashed", "dotted", or empty for solid
        };
        struct Scatter {
            std::vector<double> x, y;
            std::string color;
            std::string label;
            double radius = 3.0;
        };
        struct HLine {
            double yval;
            std::string color;
            double width = 1.0;
            std::string dash;
            std::string label;
        };
        std::vector<Line> lines;
        std::vector<Scatter> scatters;
        std::vector<HLine> hlines;
    };

    int width, height;
    int nRows = 1, nCols = 1;
    std::string title;
    std::vector<Panel> panels;

    SVGPlot(int w, int h, int rows = 1, int cols = 1)
        : width(w), height(h), nRows(rows), nCols(cols) {
        panels.resize(rows * cols);
    }

    Panel& panel(int r, int c) { return panels[r * nCols + c]; }
    Panel& panel(int i = 0) { return panels[i]; }

    void addLine(int pi, const double* x, const double* y, int n,
                 const std::string& color, const std::string& label = "",
                 double width = 1.5, const std::string& dash = "") {
        auto& p = panels[pi];
        Panel::Line l;
        l.x.assign(x, x + n); l.y.assign(y, y + n);
        l.color = color; l.label = label; l.width = width; l.dash = dash;
        p.lines.push_back(std::move(l));
        updateBounds(p);
    }

    void addLine(const double* x, const double* y, int n,
                 const std::string& color, const std::string& label = "",
                 double width = 1.5, const std::string& dash = "") {
        addLine(0, x, y, n, color, label, width, dash);
    }

    void addScatter(int pi, const double* x, const double* y, int n,
                    const std::string& color, const std::string& label = "",
                    double radius = 3.0) {
        auto& p = panels[pi];
        Panel::Scatter s;
        s.x.assign(x, x + n); s.y.assign(y, y + n);
        s.color = color; s.label = label; s.radius = radius;
        p.scatters.push_back(std::move(s));
        updateBounds(p);
    }

    void addHLine(int pi, double yval, const std::string& color,
                  double width = 1.0, const std::string& dash = "dashed",
                  const std::string& label = "") {
        auto& p = panels[pi];
        Panel::HLine h;
        h.yval = yval; h.color = color; h.width = width;
        h.dash = dash; h.label = label;
        p.hlines.push_back(h);
    }

    void setXLabel(int pi, const std::string& s) { panels[pi].xlabel = s; }
    void setYLabel(int pi, const std::string& s) { panels[pi].ylabel = s; }
    void setPanelTitle(int pi, const std::string& s) { panels[pi].title = s; }
    void setLogY(int pi, bool v = true) { panels[pi].logY = v; }

    void autoBounds(int pi) {
        auto& p = panels[pi];
        p.xmin = std::numeric_limits<double>::max();
        p.xmax = -std::numeric_limits<double>::max();
        p.ymin = p.logY ? 1e-12 : std::numeric_limits<double>::max();
        p.ymax = -std::numeric_limits<double>::max();
        for (auto& l : p.lines) {
            if (!l.x.empty()) {
                auto [a,b] = std::minmax_element(l.x.begin(), l.x.end());
                p.xmin = std::min(p.xmin, *a); p.xmax = std::max(p.xmax, *b);
            }
            if (!l.y.empty()) {
                for (double v : l.y) {
                    if (p.logY) { if (v > 0) { p.ymin = std::min(p.ymin, v); p.ymax = std::max(p.ymax, v); } }
                    else { p.ymin = std::min(p.ymin, v); p.ymax = std::max(p.ymax, v); }
                }
            }
        }
        for (auto& s : p.scatters) {
            if (!s.x.empty()) {
                auto [a,b] = std::minmax_element(s.x.begin(), s.x.end());
                p.xmin = std::min(p.xmin, *a); p.xmax = std::max(p.xmax, *b);
            }
            if (!s.y.empty()) {
                for (double v : s.y) {
                    if (p.logY) { if (v > 0) { p.ymin = std::min(p.ymin, v); p.ymax = std::max(p.ymax, v); } }
                    else { p.ymin = std::min(p.ymin, v); p.ymax = std::max(p.ymax, v); }
                }
            }
        }
        if (p.xmin >= p.xmax) { p.xmax = p.xmin + 1; }
        if (p.ymin >= p.ymax) { p.ymax = p.ymin + 1; }
        double dx = (p.xmax - p.xmin) * 0.05;
        double dy = p.logY ? 0.2 : (p.ymax - p.ymin) * 0.08;
        p.xmin -= dx; p.xmax += dx;
        if (p.logY) { p.ymin *= 0.5; p.ymax *= 2.0; }
        else { p.ymin -= dy; p.ymax += dy; }
    }

    void save(const std::string& filepath) const {
        std::ofstream f(filepath);
        if (!f.is_open()) return;

        f << R"(<?xml version="1.0" encoding="UTF-8"?>)" << "\n";
        f << "<svg xmlns='http://www.w3.org/2000/svg' "
          << "width='" << width << "' height='" << height << "' "
          << "viewBox='0 0 " << width << " " << height << "'>\n";
        f << "<rect width='100%' height='100%' fill='white'/>\n";

        // Title
        if (!title.empty()) {
            f << "<text x='" << width/2 << "' y='22' text-anchor='middle' "
              << "font-family='sans-serif' font-size='14' font-weight='bold'>"
              << escape(title) << "</text>\n";
        }

        int topMargin = title.empty() ? 10 : 32;
        int bottomMargin = 4;
        double panelW = (width - 60.0) / nCols;
        double panelH = (height - topMargin - bottomMargin - 30.0) / nRows;

        for (int pi = 0; pi < static_cast<int>(panels.size()); ++pi) {
            int r = pi / nCols;
            int c = pi % nCols;

            double px = 55 + c * panelW;
            double py = topMargin + r * panelH;
            double pw = panelW - 15;
            double ph = panelH - 30;

            renderPanel(f, panels[pi], px, py, pw, ph);
        }

        f << "</svg>\n";
    }

private:
    void updateBounds(Panel& p) const {
        // bounds set manually or via autoBounds()
    }

    static std::string escape(const std::string& s) {
        std::string o;
        for (char ch : s) {
            if (ch == '<') o += "&lt;";
            else if (ch == '>') o += "&gt;";
            else if (ch == '&') o += "&amp;";
            else o += ch;
        }
        return o;
    }

    // SVG y is top-down; our data y is bottom-up.
    double mapX(double x, double xmin, double xmax, double px, double pw) const {
        return px + pw * (x - xmin) / (xmax - xmin);
    }
    double mapY(double y, double ymin, double ymax, double py, double ph, bool logY) const {
        if (logY) {
            if (y <= 0) y = ymin;
            return py + ph - ph * (std::log10(y) - std::log10(ymin))
                     / (std::log10(ymax) - std::log10(ymin));
        }
        return py + ph - ph * (y - ymin) / (ymax - ymin);
    }

    void renderPanel(std::ofstream& f, const Panel& p,
                     double px, double py, double pw, double ph) const {
        double xmin = p.xmin, xmax = p.xmax;
        double ymin = p.ymin, ymax = p.ymax;

        f << "<g>\n";

        // Panel background
        f << "<rect x='" << px << "' y='" << py << "' width='" << pw
          << "' height='" << ph << "' fill='#fafafa' stroke='#ccc' stroke-width='1'/>\n";

        // Grid lines
        int nGrid = 5;
        for (int i = 0; i <= nGrid; ++i) {
            double t = static_cast<double>(i) / nGrid;
            double gx = px + t * pw;
            double gy = py + ph - t * ph;
            f << "<line x1='" << gx << "' y1='" << py << "' x2='" << gx
              << "' y2='" << py+ph << "' stroke='#e0e0e0' stroke-width='0.5'/>\n";
            f << "<line x1='" << px << "' y1='" << gy << "' x2='" << px+pw
              << "' y2='" << gy << "' stroke='#e0e0e0' stroke-width='0.5'/>\n";
        }

        // Clip region for data
        std::string clipId = "clip_" + std::to_string(reinterpret_cast<uintptr_t>(&p));
        f << "<defs><clipPath id='" << clipId << "'>"
          << "<rect x='" << px << "' y='" << py << "' width='" << pw
          << "' height='" << ph << "'/>"
          << "</clipPath></defs>\n";
        f << "<g clip-path='url(#" << clipId << ")'>\n";

        // Horizontal lines
        for (auto& hl : p.hlines) {
            double hy = mapY(hl.yval, ymin, ymax, py, ph, p.logY);
            std::string dashAttr;
            if (hl.dash == "dashed") dashAttr = " stroke-dasharray='6,4'";
            else if (hl.dash == "dotted") dashAttr = " stroke-dasharray='2,3'";
            f << "<line x1='" << px << "' y1='" << hy << "' x2='" << px+pw
              << "' y2='" << hy << "' stroke='" << hl.color
              << "' stroke-width='" << hl.width << "'" << dashAttr << "/>\n";
            if (!hl.label.empty()) {
                f << "<text x='" << (px+pw-4) << "' y='" << (hy-4)
                  << "' fill='" << hl.color << "' font-family='sans-serif' font-size='9' "
                  << "text-anchor='end'>" << escape(hl.label) << "</text>\n";
            }
        }

        // Data lines
        for (auto& l : p.lines) {
            if (l.x.size() < 2) continue;
            std::string dashAttr;
            if (l.dash == "dashed") dashAttr = " stroke-dasharray='6,4'";
            else if (l.dash == "dotted") dashAttr = " stroke-dasharray='2,3'";

            f << "<polyline fill='none' stroke='" << l.color
              << "' stroke-width='" << l.width << "'"
              << dashAttr << " points='";
            for (size_t i = 0; i < l.x.size(); ++i) {
                double sx = mapX(l.x[i], xmin, xmax, px, pw);
                double sy = mapY(l.y[i], ymin, ymax, py, ph, p.logY);
                f << sx << "," << sy << " ";
            }
            f << "'/>\n";
        }

        // Scatter points
        for (auto& s : p.scatters) {
            for (size_t i = 0; i < s.x.size(); ++i) {
                double sx = mapX(s.x[i], xmin, xmax, px, pw);
                double sy = mapY(s.y[i], ymin, ymax, py, ph, p.logY);
                f << "<circle cx='" << sx << "' cy='" << sy
                  << "' r='" << s.radius << "' fill='" << s.color << "'/>\n";
            }
        }

        f << "</g>\n";  // end clip

        // Axis lines
        double axY = mapY(ymin, ymin, ymax, py, ph, p.logY);
        f << "<line x1='" << px << "' y1='" << axY << "' x2='" << px+pw
          << "' y2='" << axY << "' stroke='#333' stroke-width='1'/>\n";
        f << "<line x1='" << px << "' y1='" << py << "' x2='" << px
          << "' y2='" << py+ph << "' stroke='#333' stroke-width='1'/>\n";

        // Tick marks & labels on Y axis
        for (int i = 0; i <= 4; ++i) {
            double t = i / 4.0;
            double val = p.logY
                ? ymin * std::pow(ymax/ymin, t)
                : ymin + t * (ymax - ymin);
            double sy = mapY(val, ymin, ymax, py, ph, p.logY);
            f << "<line x1='" << (px-4) << "' y1='" << sy << "' x2='" << px
              << "' y2='" << sy << "' stroke='#333' stroke-width='1'/>\n";

            std::stringstream ss;
            if (p.logY) ss << std::scientific << std::setprecision(1) << val;
            else if (std::abs(val) < 0.01 || std::abs(val) >= 1e4)
                ss << std::scientific << std::setprecision(1) << val;
            else ss << std::fixed << std::setprecision(2) << val;

            f << "<text x='" << (px-7) << "' y='" << (sy+4)
              << "' fill='#333' font-family='sans-serif' font-size='9' "
              << "text-anchor='end'>" << escape(ss.str()) << "</text>\n";
        }

        // Tick marks & labels on X axis
        for (int i = 0; i <= 4; ++i) {
            double t = i / 4.0;
            double val = xmin + t * (xmax - xmin);
            double sx = mapX(val, xmin, xmax, px, pw);
            f << "<line x1='" << sx << "' y1='" << axY << "' x2='" << sx
              << "' y2='" << (axY+4) << "' stroke='#333' stroke-width='1'/>\n";

            std::stringstream ss;
            ss << std::fixed << std::setprecision(0) << val;
            f << "<text x='" << sx << "' y='" << (axY+15)
              << "' fill='#333' font-family='sans-serif' font-size='9' "
              << "text-anchor='middle'>" << escape(ss.str()) << "</text>\n";
        }

        // Axis labels
        if (!p.xlabel.empty()) {
            f << "<text x='" << (px+pw/2) << "' y='" << (py+ph+28)
              << "' fill='#333' font-family='sans-serif' font-size='11' "
              << "text-anchor='middle'>" << escape(p.xlabel) << "</text>\n";
        }
        if (!p.ylabel.empty()) {
            f << "<text transform='translate(" << (px-42) << "," << (py+ph/2)
              << ") rotate(-90)' fill='#333' font-family='sans-serif' "
              << "font-size='11' text-anchor='middle'>"
              << escape(p.ylabel) << "</text>\n";
        }

        // Panel title
        if (!p.title.empty()) {
            f << "<text x='" << (px+pw/2) << "' y='" << (py-5)
              << "' fill='#222' font-family='sans-serif' font-size='11' "
              << "font-weight='bold' text-anchor='middle'>"
              << escape(p.title) << "</text>\n";
        }

        // Legend
        bool haveLegend = false;
        for (auto& l : p.lines) if (!l.label.empty()) { haveLegend = true; break; }
        for (auto& s : p.scatters) if (!s.label.empty()) { haveLegend = true; break; }
        for (auto& hl : p.hlines) if (!hl.label.empty()) { haveLegend = true; break; }

        if (haveLegend) {
            double lx = px + pw - 160, ly = py + 6;
            f << "<rect x='" << lx << "' y='" << ly << "' width='155' height='"
              << (18 * (countLabels(p)+1)) << "' fill='white' fill-opacity='0.85' "
              << "stroke='#ddd' stroke-width='0.5' rx='3'/>\n";

            int li = 0;
            for (auto& l : p.lines) {
                if (l.label.empty()) continue;
                double lly = ly + 16 + li * 18;
                f << "<line x1='" << (lx+8) << "' y1='" << lly << "' x2='"
                  << (lx+38) << "' y2='" << lly << "' stroke='" << l.color
                  << "' stroke-width='2'/>\n";
                f << "<text x='" << (lx+44) << "' y='" << (lly+4)
                  << "' fill='#333' font-family='sans-serif' font-size='9'>"
                  << escape(l.label) << "</text>\n";
                ++li;
            }
            for (auto& s : p.scatters) {
                if (s.label.empty()) continue;
                double lly = ly + 16 + li * 18;
                f << "<circle cx='" << (lx+23) << "' cy='" << lly
                  << "' r='3.5' fill='" << s.color << "'/>\n";
                f << "<text x='" << (lx+44) << "' y='" << (lly+4)
                  << "' fill='#333' font-family='sans-serif' font-size='9'>"
                  << escape(s.label) << "</text>\n";
                ++li;
            }
            for (auto& hl : p.hlines) {
                if (hl.label.empty()) continue;
                double lly = ly + 16 + li * 18;
                f << "<line x1='" << (lx+8) << "' y1='" << lly << "' x2='"
                  << (lx+38) << "' y2='" << lly << "' stroke='" << hl.color
                  << "' stroke-width='1.5' stroke-dasharray='6,4'/>\n";
                f << "<text x='" << (lx+44) << "' y='" << (lly+4)
                  << "' fill='#333' font-family='sans-serif' font-size='9'>"
                  << escape(hl.label) << "</text>\n";
                ++li;
            }
        }

        f << "</g>\n";
    }

    static int countLabels(const Panel& p) {
        int n = 0;
        for (auto& l : p.lines) if (!l.label.empty()) ++n;
        for (auto& s : p.scatters) if (!s.label.empty()) ++n;
        for (auto& h : p.hlines) if (!h.label.empty()) ++n;
        return n;
    }
};

} // namespace litho_invert
