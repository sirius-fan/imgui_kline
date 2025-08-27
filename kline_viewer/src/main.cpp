// Use ImGui OpenGL3 loader to avoid extra deps
#include "imgui_impl_opengl3_loader.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>

#include <GLFW/glfw3.h>

#include "indicators.hpp"

struct ViewState {
    float scale_x = 6.0f;     // pixels per candle
    float scale_y = 0.5f;     // pixels per price unit (dynamic)
    float scroll_x = 0.0f;    // index offset
    float top_padding = 16.0f;
    float bottom_padding = 80.0f; // room for MACD/RSI
};

struct ChartOptions {
    bool show_sma20 = true;
    bool show_ema50 = true;
    bool show_macd = true;
    bool show_rsi = true;
};

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

static std::vector<Candle> gen_sample_data(size_t n = 600) {
    std::vector<Candle> v; v.reserve(n);
    double price = 100.0;
    auto t0 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double drift = (std::sin(i*0.03) + std::sin(i*0.013)) * 0.5;
        double vol = 0.8 + 0.4*std::sin(i*0.17 + 0.5);
        double open = price;
        double high = open + std::abs(drift*2.0 + 0.5)*vol;
        double low  = open - std::abs(drift*2.0 + 0.5)*vol;
        double close = open + drift * vol;
        double volume = 1000.0 + 400.0 * std::sin(i*0.07);
        price = close;
        v.push_back({t0 + (double)i, open, high, low, close, volume});
    }
    return v;
}

static void draw_grid(const ImVec2& p0, const ImVec2& p1, float y_min, float y_max) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImU32 col_grid = IM_COL32(64,64,64,80);
    // Horizontal lines
    int h_lines = 6;
    for (int i=0; i<=h_lines; ++i) {
        float t = (float)i / (float)h_lines;
    float y = p0.y + (p1.y - p0.y) * t;
        dl->AddLine(ImVec2(p0.x, y), ImVec2(p1.x, y), col_grid);
    float price = y_max + (y_min - y_max) * t;
        char buf[32];
        snprintf(buf, sizeof(buf), "%.2f", price);
        dl->AddText(ImVec2(p0.x + 4, y - ImGui::GetTextLineHeight()*0.5f), IM_COL32(170,170,170,200), buf);
    }
}

static void draw_candles(const std::vector<Candle>& data, const ViewState& vs, const ImVec2& canvas_pos, const ImVec2& canvas_size,
                         int begin, int end, float y_min, float y_max) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float W = canvas_size.x;
    const float H = canvas_size.y;
    const float candle_w = std::max(1.0f, vs.scale_x * 0.7f);

    for (int i=begin; i<end; ++i) {
        const Candle& c = data[i];
        float x = canvas_pos.x + (i - vs.scroll_x) * vs.scale_x;
        float x0 = x - candle_w * 0.5f;
        float x1 = x + candle_w * 0.5f;
        auto y_to_screen = [&](double y){
            float ty = (float)((y - y_min) / (y_max - y_min));
            return canvas_pos.y + (1.0f - ty) * H;
        };
        float y_open = y_to_screen(c.open);
        float y_close = y_to_screen(c.close);
        float y_high = y_to_screen(c.high);
        float y_low  = y_to_screen(c.low);
        bool up = c.close >= c.open;
        ImU32 col = up ? IM_COL32(82, 196, 26, 255) : IM_COL32(255, 77, 79, 255);
        // Wick
        dl->AddLine(ImVec2(x, y_high), ImVec2(x, y_low), col, 1.0f);
        // Body
        if (std::abs(y_open - y_close) < 1.0f)
            dl->AddLine(ImVec2(x0, y_open), ImVec2(x1, y_close), col, candle_w);
        else
            dl->AddRectFilled(ImVec2(x0, std::min(y_open,y_close)), ImVec2(x1, std::max(y_open,y_close)), col);
    }
}

static void draw_line_series(const std::vector<double>& series, const ViewState& vs,
                             const ImVec2& canvas_pos, const ImVec2& canvas_size,
                             int begin, int end, float y_min, float y_max, ImU32 col) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float H = canvas_size.y;
    auto y_to_screen = [&](double y){
        float ty = (float)((y - y_min) / (y_max - y_min));
        return canvas_pos.y + (1.0f - ty) * H;
    };
    ImVec2 prev;
    bool has_prev = false;
    for (int i=begin; i<end; ++i) {
        double v = series[i];
        if (std::isnan(v)) { has_prev = false; continue; }
        float x = canvas_pos.x + (i - vs.scroll_x) * vs.scale_x;
        float y = y_to_screen(v);
        ImVec2 cur(x, y);
        if (has_prev) dl->AddLine(prev, cur, col, 1.5f);
        prev = cur; has_prev = true;
    }
}

int main(int, char**)
{
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    GLFWwindow* window = glfwCreateWindow(1280, 800, "ImGui K-Line Viewer", NULL, NULL);
    if (window == NULL)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    // Docking is not available on master branch by default
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Data
    std::vector<Candle> candles = gen_sample_data(800);
    std::vector<double> closes; closes.reserve(candles.size());
    for (auto& c: candles) closes.push_back(c.close);

    std::vector<double> sma20 = ind::sma(closes, 20);
    std::vector<double> ema50 = ind::ema(closes, 50);
    std::vector<double> macd_line, signal_line, hist;
    ind::macd(closes, 12, 26, 9, &macd_line, &signal_line, &hist);
    std::vector<double> rsi14; ind::rsi(closes, 14, rsi14);

    ViewState vs;
    ChartOptions opt;
    bool crosshair_visible = false;
    int cross_idx = -1;
    double cross_price = 0.0;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

    // No docking on master branch

        ImGui::Begin("Chart");
    ImVec2 canvas_pos = ImGui::GetCursorScreenPos();
    ImVec2 canvas_size = ImGui::GetContentRegionAvail();
    if (canvas_size.x < 50) canvas_size.x = 50;
    if (canvas_size.y < 200) canvas_size.y = 200; // need room for subplots

    // Layout: main, optional RSI, optional MACD stacked vertically
    const float spacing = 6.0f;
    const float rsi_h = opt.show_rsi ? 80.0f : 0.0f;
    const float macd_h = opt.show_macd ? 120.0f : 0.0f;
    float main_h = canvas_size.y - (opt.show_rsi ? (rsi_h + spacing) : 0.0f)
                      - (opt.show_macd ? (macd_h + spacing) : 0.0f);
    if (main_h < 120.0f) main_h = 120.0f;
    ImVec2 main_pos = canvas_pos;
    ImVec2 main_size = ImVec2(canvas_size.x, main_h);
    ImVec2 rsi_pos = ImVec2(main_pos.x, main_pos.y + main_size.y + spacing);
    ImVec2 rsi_size = ImVec2(canvas_size.x, rsi_h);
    ImVec2 macd_pos = opt.show_rsi ? ImVec2(rsi_pos.x, rsi_pos.y + rsi_size.y + spacing)
                       : ImVec2(main_pos.x, main_pos.y + main_size.y + spacing);
    ImVec2 macd_size = ImVec2(canvas_size.x, macd_h);

        // Interaction: mouse wheel zoom/scroll
        ImGuiIO& io = ImGui::GetIO();
        if (ImGui::IsWindowHovered()) {
            float wheel = io.MouseWheel;
            if (wheel != 0.0f) {
                float mouse_x = io.MousePos.x - main_pos.x;
                float center_index = vs.scroll_x + mouse_x / vs.scale_x;
                float prev_scale = vs.scale_x;
                vs.scale_x = std::clamp(vs.scale_x * (1.0f + wheel*0.1f), 1.5f, 30.0f);
                vs.scroll_x = center_index - mouse_x / vs.scale_x;
            }
            // right click drag
            if (ImGui::IsMouseDown(1) && ImGui::IsMouseDragging(1)) {
                float dx = io.MouseDelta.x;
                vs.scroll_x -= dx / vs.scale_x;
                vs.scroll_x = std::clamp(vs.scroll_x, 0.0f, (float)candles.size());
            }
            // left click toggles crosshair
            if (ImGui::IsMouseClicked(0)) {
                // Only toggle when clicking inside the chart region
                ImVec2 mp = io.MousePos;
                bool in_main = mp.x >= main_pos.x && mp.x <= main_pos.x + main_size.x && mp.y >= main_pos.y && mp.y <= main_pos.y + main_size.y;
                bool in_rsi  = opt.show_rsi && mp.x >= rsi_pos.x && mp.x <= rsi_pos.x + rsi_size.x && mp.y >= rsi_pos.y && mp.y <= rsi_pos.y + rsi_size.y;
                bool in_macd = opt.show_macd && mp.x >= macd_pos.x && mp.x <= macd_pos.x + macd_size.x && mp.y >= macd_pos.y && mp.y <= macd_pos.y + macd_size.y;
                if (in_main || in_rsi || in_macd) crosshair_visible = !crosshair_visible;
            }
        }

        // Determine visible range
        int begin = std::max(0, (int)std::floor(vs.scroll_x));
        int count = (int)std::ceil(canvas_size.x / vs.scale_x) + 2;
        int end = std::min((int)candles.size(), begin + count);

        // Compute y-range for visible candles
        double y_min = 1e9, y_max = -1e9;
        for (int i=begin; i<end; ++i) {
            y_min = std::min(y_min, candles[i].low);
            y_max = std::max(y_max, candles[i].high);
        }
        float pad = (float)((y_max - y_min) * 0.05 + 1e-6);
        y_min -= pad; y_max += pad;

    // Draw backgrounds and borders per panel
    ImDrawList* dl = ImGui::GetWindowDrawList();
    auto rect = [&](ImVec2 p, ImVec2 s){ dl->AddRectFilled(p, ImVec2(p.x+s.x, p.y+s.y), IM_COL32(20,20,20,255)); dl->AddRect(p, ImVec2(p.x+s.x, p.y+s.y), IM_COL32(100,100,100,255)); };
    rect(main_pos, main_size);
    if (opt.show_rsi)  dl->AddRect(rsi_pos,  ImVec2(rsi_pos.x+rsi_size.x,   rsi_pos.y+rsi_size.y),   IM_COL32(100,100,100,150));
    if (opt.show_macd) dl->AddRect(macd_pos, ImVec2(macd_pos.x+macd_size.x, macd_pos.y+macd_size.y), IM_COL32(100,100,100,150));

    // Main chart
    draw_grid(main_pos, ImVec2(main_pos.x+main_size.x, main_pos.y+main_size.y), (float)y_min, (float)y_max);
    draw_candles(candles, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max);
    if (opt.show_sma20) draw_line_series(sma20, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, IM_COL32(255, 193, 7, 255));
    if (opt.show_ema50) draw_line_series(ema50, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, IM_COL32(24, 144, 255, 255));

        // MACD panel
        if (opt.show_macd) {
            // Scale MACD
            double mmin=1e9, mmax=-1e9; for (int i=begin;i<end;++i){ mmin=std::min({mmin, macd_line[i], signal_line[i], hist[i]}); mmax=std::max({mmax, macd_line[i], signal_line[i], hist[i]}); }
            mmin = std::min(mmin, 0.0); mmax = std::max(mmax, 0.0);
            auto macd_y = [&](double v){ float t=(float)((v - mmin) / (mmax - mmin + 1e-9)); return macd_pos.y + (1.0f - t) * macd_size.y; };
            // Zero line
            float y0 = macd_y(0.0);
            dl->AddLine(ImVec2(macd_pos.x, y0), ImVec2(macd_pos.x+macd_size.x, y0), IM_COL32(120,120,120,180));
            // Histogram
            for (int i=begin;i<end;++i) {
                if (std::isnan(hist[i])) continue;
                float x = macd_pos.x + (i - vs.scroll_x) * vs.scale_x;
                float x0 = x - std::max(1.0f, vs.scale_x*0.6f)/2.0f;
                float x1 = x + std::max(1.0f, vs.scale_x*0.6f)/2.0f;
                float yv = macd_y(hist[i]);
                ImU32 col = hist[i] >= 0 ? IM_COL32(0,200,0,180) : IM_COL32(220,0,0,180);
                dl->AddRectFilled(ImVec2(x0, y0), ImVec2(x1, yv), col);
            }
            // Lines
            auto draw_line_local = [&](const std::vector<double>& s, ImU32 col){
                ImVec2 p0 = macd_pos, sz = macd_size;
                const float H = sz.y;
                auto y_to = [&](double v){ float t=(float)((v - mmin) / (mmax - mmin + 1e-9)); return p0.y + (1.0f - t) * H; };
                ImVec2 prev; bool hp=false;
                for (int i=begin;i<end;++i){ double vv=s[i]; if (std::isnan(vv)) {hp=false; continue;} float x=p0.x + (i - vs.scroll_x)*vs.scale_x; float y=y_to(vv); ImVec2 cur(x,y); if(hp) dl->AddLine(prev, cur, col, 1.5f); prev=cur; hp=true; }
            };
            draw_line_local(macd_line, IM_COL32(255,255,255,220));
            draw_line_local(signal_line, IM_COL32(255,215,0,220));
        }

        // RSI panel
        if (opt.show_rsi) {
            auto rsi_y = [&](double v){ return rsi_pos.y + (float)( (100.0 - v) / 100.0 ) * rsi_size.y; };
            // 30/70 bands
            dl->AddLine(ImVec2(rsi_pos.x, rsi_y(30)), ImVec2(rsi_pos.x+rsi_size.x, rsi_y(30)), IM_COL32(150,150,150,180));
            dl->AddLine(ImVec2(rsi_pos.x, rsi_y(70)), ImVec2(rsi_pos.x+rsi_size.x, rsi_y(70)), IM_COL32(150,150,150,180));
            // RSI line
            ImVec2 prev; bool has_prev=false;
            for (int i=begin;i<end;++i){ double v=rsi14[i]; if (std::isnan(v)) {has_prev=false; continue;} float x=rsi_pos.x + (i - vs.scroll_x)*vs.scale_x; float y=rsi_y(v); ImVec2 cur(x,y); if(has_prev) dl->AddLine(prev, cur, IM_COL32(64,158,255,255), 1.5f); prev=cur; has_prev=true; }
        }

        // Crosshair on main panel and data readout
        auto format_opt = [](double v, char* buf, size_t n){ if (std::isnan(v)) { snprintf(buf,n,"-"); } else { snprintf(buf,n,"%.4f", v); } };
        if (crosshair_visible) {
            ImVec2 mp = io.MousePos;
            // Snap to nearest candle index based on mouse x
            float mouse_x_rel = std::clamp(mp.x - main_pos.x, 0.0f, main_size.x);
            int idx = (int)std::round(vs.scroll_x + mouse_x_rel / vs.scale_x);
            idx = std::clamp(idx, 0, (int)candles.size()-1);
            cross_idx = idx;
            // Compute price at mouse y
            float mouse_y_rel = std::clamp(mp.y - main_pos.y, 0.0f, main_size.y);
            cross_price = (double) (y_max - (mouse_y_rel / main_size.y) * (y_max - y_min));

            // Crosshair lines
            float cx = main_pos.x + (cross_idx - vs.scroll_x) * vs.scale_x;
            float cy_main = main_pos.y + (float)((y_max - cross_price) / (y_max - y_min)) * main_size.y;
            dl->AddLine(ImVec2(cx, main_pos.y), ImVec2(cx, main_pos.y + main_size.y), IM_COL32(200,200,200,120));
            dl->AddLine(ImVec2(main_pos.x, cy_main), ImVec2(main_pos.x + main_size.x, cy_main), IM_COL32(200,200,200,120));
            if (opt.show_rsi)  dl->AddLine(ImVec2(cx, rsi_pos.y),  ImVec2(cx, rsi_pos.y  + rsi_size.y),  IM_COL32(200,200,200,60));
            if (opt.show_macd) dl->AddLine(ImVec2(cx, macd_pos.y), ImVec2(cx, macd_pos.y + macd_size.y), IM_COL32(200,200,200,60));

            // Data box
            const Candle& c = candles[cross_idx];
            char buf1[256];
            char sv[32], ev[32], mv[32], sg[32], hs[32], rs[32];
            format_opt(sma20[cross_idx], sv, sizeof(sv));
            format_opt(ema50[cross_idx], ev, sizeof(ev));
            format_opt(macd_line[cross_idx], mv, sizeof(mv));
            format_opt(signal_line[cross_idx], sg, sizeof(sg));
            format_opt(hist[cross_idx], hs, sizeof(hs));
            format_opt(rsi14[cross_idx], rs, sizeof(rs));
            snprintf(buf1, sizeof(buf1),
                     "Idx: %d\nO: %.4f  H: %.4f  L: %.4f  C: %.4f\nSMA20: %s  EMA50: %s\nMACD: %s  Signal: %s  Hist: %s\nRSI14: %s",
                     cross_idx, c.open, c.high, c.low, c.close,
                     sv, ev, mv, sg, hs, rs);
            ImVec2 box_pos = ImVec2(main_pos.x + 8, main_pos.y + 8);
            ImGui::SetCursorScreenPos(box_pos);
            ImGui::BeginChild("DataBox", ImVec2(320, ImGui::GetTextLineHeightWithSpacing()*5.5f), ImGuiChildFlags_Border | ImGuiChildFlags_AutoResizeY);
            ImGui::TextUnformatted(buf1);
            ImGui::EndChild();
        }

        // UI controls
        ImGui::SetCursorScreenPos(ImVec2(canvas_pos.x+8, canvas_pos.y+8));
    ImGui::BeginChild("Legend", ImVec2(260,160), ImGuiChildFlags_Border);
        ImGui::Text("K-Line Viewer");
        ImGui::Text("Candles: %zu", candles.size());
        ImGui::SliderFloat("Scale X", &vs.scale_x, 1.5f, 30.0f);
        ImGui::Text("Scroll: %.1f", vs.scroll_x);
    ImGui::Separator();
    ImGui::Text("Indicators");
    ImGui::Checkbox("SMA20", &opt.show_sma20);
    ImGui::Checkbox("EMA50", &opt.show_ema50);
    ImGui::Checkbox("MACD", &opt.show_macd);
    ImGui::Checkbox("RSI", &opt.show_rsi);
    ImGui::Separator();
    ImGui::Text("Crosshair: %s (L-Click to toggle)", crosshair_visible ? "ON" : "OFF");
        ImGui::EndChild();

        ImGui::End();

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
