// Use ImGui OpenGL3 loader to avoid extra deps
#include "imgui_impl_opengl3_loader.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <cstdint>
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstring>
#include <limits>

#include <GLFW/glfw3.h>
// SQLite for DB data source
#include <sqlite3.h>

#include "indicators.hpp"

// Find nearest candle index by time using binary search (candles sorted by time ascending)
static int find_nearest_index_by_time(const std::vector<Candle> &cs, std::time_t t) {
    if (cs.empty()) return -1;
    int lo = 0, hi = (int)cs.size();
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        std::time_t mt = (std::time_t)cs[mid].time;
        if (mt < t)
            lo = mid + 1;
        else
            hi = mid;
    }
    int pos = lo;
    if (pos >= (int)cs.size()) return (int)cs.size() - 1;
    if (pos > 0) {
        std::time_t t1 = (std::time_t)cs[pos].time;
        std::time_t t0 = (std::time_t)cs[pos - 1].time;
        long long d0 = std::llabs((long long)t - (long long)t0);
        long long d1 = std::llabs((long long)t - (long long)t1);
        return (d0 <= d1) ? (pos - 1) : pos;
    }
    return pos;
}

struct ViewState {
    float scale_x = 6.0f;  // pixels per candle
    float scale_y = 0.5f;  // pixels per price unit (dynamic)
    float scroll_x = 0.0f; // index offset
    float top_padding = 16.0f;
    float bottom_padding = 80.0f; // room for MACD/RSI
};

struct ChartOptions {
    bool show_sma20 = true;
    bool show_ema50 = true;
    bool show_macd = true;
    bool show_rsi = true;
    bool show_volume = true;
    bool labels_follow_cursor = false; // when true, right-side floating labels follow crosshair index
    bool show_close_line = false;      // show close price as a line on main chart
    bool show_boll = false;            // show Bollinger Bands
    bool show_hlc_area = false;        // show HLC area (High/Close/Low with fills)
    bool show_kdj = false;             // show KDJ subpanel
    bool show_sar = false;             // show Parabolic SAR on main chart
    bool show_td9 = false;             // show TD Sequential Setup (1..9) on main chart
    // BOLL styles
    ImVec4 boll_mid_color = ImVec4(0.77f, 0.35f, 0.94f, 0.86f); // ~ IM_COL32(197,90,240,220)
    ImVec4 boll_band_color = ImVec4(0.77f, 0.35f, 0.94f, 0.55f); // ~ IM_COL32(197,90,240,140)
    float  boll_thickness = 1.6f;
    // SAR style
    float sar_dot_size = 3.0f;
};

// Trade annotation for buy/sell pair
struct TradeAnnot {
    std::time_t buy_t{};
    double buy_p{};
    std::time_t sell_t{};
    double sell_p{};
    bool connect_line{true};
};

static void glfw_error_callback(int error, const char *description) {
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

static std::vector<Candle> gen_sample_data(size_t n = 600) {
    std::vector<Candle> v;
    v.reserve(n);
    double price = 100.0;
    auto t0 = 0;
    for (size_t i = 0; i < n; ++i) {
        double drift = (std::sin(i * 0.03) + std::sin(i * 0.013)) * 0.5;
        double vol = 0.8 + 0.4 * std::sin(i * 0.17 + 0.5);
        double open = price;
        double high = open + std::abs(drift * 2.0 + 0.5) * vol;
        double low = open - std::abs(drift * 2.0 + 0.5) * vol;
        double close = open + drift * vol;
        double volume = 1000.0 + 400.0 * std::sin(i * 0.07);
        price = close;
        v.push_back({t0 + (uint64_t)i, open, high, low, close, volume});
    }
    return v;
}

static inline bool parse_datetime_to_time_t(const std::string &date, const std::string &time, std::time_t &out) {
    // date: YYYY-MM-DD, time: HH:MM (24h)
    int Y = 0, M = 0, D = 0, h = 0, m = 0;
    if (sscanf(date.c_str(), "%d-%d-%d", &Y, &M, &D) != 3) return false;
    if (sscanf(time.c_str(), "%d:%d", &h, &m) != 2) {
        h = 0;
        m = 0;
    }
    std::tm tm{};
    tm.tm_year = Y - 1900;
    tm.tm_mon = M - 1;
    tm.tm_mday = D;
    tm.tm_hour = h;
    tm.tm_min = m;
    tm.tm_sec = 0;
    // Use local time; for pure ordering this is fine
    out = std::mktime(&tm);
    return out != (std::time_t)-1;
}

static bool load_csv_sh_index(const std::string &path, std::vector<Candle> &out) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) return false;
    std::string line;
    // Skip header
    if (!std::getline(ifs, line)) return false;
    out.clear();
    out.reserve(4096);
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::vector<std::string> cols;
        cols.reserve(12);
        std::string cur;
        cur.reserve(32);
        std::istringstream ss(line);
        while (std::getline(ss, cur, ',')) cols.push_back(cur);
        if (cols.size() < 10) continue; // 日期,时间,开盘,最高,最低,收盘,成交量,成交额,涨跌价,涨跌幅
        std::time_t tt{};
        if (!parse_datetime_to_time_t(cols[0], cols[1], tt)) continue;
        auto to_d = [](const std::string &s) { if (s.empty()) return 0.0; try { return std::stod(s); } catch(...) { return 0.0; } };
        double open = to_d(cols[2]);
        double high = to_d(cols[3]);
        double low = to_d(cols[4]);
        double close = to_d(cols[5]);
        double vol = to_d(cols[6]);
        out.push_back({(uint64_t)tt, open, high, low, close, vol});
    }
    return !out.empty();
}

static inline void format_time_label(uint64_t t, char *buf, size_t n, bool with_time = false) {
    std::time_t tt = (std::time_t)t;
    std::tm *lt = std::localtime(&tt);
    if (!lt) {
        snprintf(buf, n, "-");
        return;
    }
    if (with_time)
        std::strftime(buf, n, "%Y-%m-%d %H:%M", lt);
    else
        std::strftime(buf, n, "%Y-%m-%d", lt);
}

// --- SQLite 1m loader ---
struct DbLoadOptions {
    std::string db_path;
    std::string symbol;
    int64_t start_ms = 0;             // inclusive
    int64_t end_ms = (int64_t)9e18;   // inclusive
};

static bool load_sqlite_1m(const DbLoadOptions &opt, std::vector<Candle> &out) {
    out.clear();
    sqlite3 *db = nullptr;
    int rc = sqlite3_open(opt.db_path.c_str(), &db);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQLite open failed: %s\n", sqlite3_errmsg(db));
        if (db) sqlite3_close(db);
        return false;
    }
    const char *sql = "SELECT open_time, symbol, close_time, open_price, high_price, low_price, close_price, volume "
                      "FROM klines_1m WHERE symbol = ? AND open_time >= ? AND open_time <= ? ORDER BY open_time ASC;";
    sqlite3_stmt *stmt = nullptr;
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQLite prepare failed: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return false;
    }
    sqlite3_bind_text(stmt, 1, opt.symbol.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int64(stmt, 2, opt.start_ms);
    sqlite3_bind_int64(stmt, 3, opt.end_ms);
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        int64_t open_time_ms = sqlite3_column_int64(stmt, 0);
        // const unsigned char* sym = sqlite3_column_text(stmt, 1);
        // int64_t close_time_ms = sqlite3_column_int64(stmt, 2);
        double open = sqlite3_column_double(stmt, 3);
        double high = sqlite3_column_double(stmt, 4);
        double low  = sqlite3_column_double(stmt, 5);
        double close = sqlite3_column_double(stmt, 6);
        double vol = sqlite3_column_double(stmt, 7);
        Candle c;
        c.time = (uint64_t)(open_time_ms / 1000); // seconds
        c.open = open; c.high = high; c.low = low; c.close = close; c.volume = vol;
        out.push_back(c);
    }
    bool ok = (rc == SQLITE_DONE);
    if (!ok) fprintf(stderr, "SQLite step error: %s\n", sqlite3_errmsg(db));
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    return ok && !out.empty();
}

// --- Aggregate 1m to higher timeframe (O/H/L/C/V) ---
static std::vector<Candle> aggregate_timeframe(const std::vector<Candle> &base_1m, int interval_sec) {
    std::vector<Candle> out;
    if (base_1m.empty() || interval_sec <= 60) return base_1m; // 1m passthrough
    Candle cur{};
    bool has = false;
    uint64_t bucket = 0;
    for (const auto &c : base_1m) {
        uint64_t b = (c.time / (uint64_t)interval_sec) * (uint64_t)interval_sec;
        if (!has) {
            has = true;
            bucket = b;
            cur.time = bucket;
            cur.open = c.open;
            cur.high = c.high;
            cur.low  = c.low;
            cur.close = c.close;
            cur.volume = c.volume;
        } else if (b == bucket) {
            cur.high = std::max(cur.high, c.high);
            cur.low  = std::min(cur.low, c.low);
            cur.close = c.close;
            cur.volume += c.volume;
        } else {
            out.push_back(cur);
            bucket = b;
            cur.time = bucket;
            cur.open = c.open;
            cur.high = c.high;
            cur.low  = c.low;
            cur.close = c.close;
            cur.volume = c.volume;
        }
    }
    if (has) out.push_back(cur);
    return out;
}

/**
 * 绘制面板内的价格网格与刻度文本（仅水平线）。
 *
 * 参数说明：
 * - p0: 面板左上角坐标（屏幕坐标，像素）
 * - p1: 面板右下角坐标（屏幕坐标，像素）
 * - y_min/y_max: 当前面板对应的数据值范围，y 轴从上到下映射 y_max → y_min
 *
 * 说明：
 * - 仅负责网格与刻度文本的绘制，不包含裁剪；调用方应在调用前设置好裁剪矩形。
 * - 刻度文本使用 y_max→y_min 的线性映射计算对应的价格值。
 */
static void draw_grid(const ImVec2 &p0, const ImVec2 &p1, float y_min, float y_max) {
    ImDrawList *dl = ImGui::GetWindowDrawList();
    const ImU32 col_grid = IM_COL32(64, 64, 64, 80);
    // Horizontal lines
    int h_lines = 6;
    for (int i = 0; i <= h_lines; ++i) {
        float t = (float)i / (float)h_lines;
        float y = p0.y + (p1.y - p0.y) * t;
        dl->AddLine(ImVec2(p0.x, y), ImVec2(p1.x, y), col_grid);
        float price = y_max + (y_min - y_max) * t;
        char buf[32];
        snprintf(buf, sizeof(buf), "%.2f", price);
        dl->AddText(ImVec2(p0.x + 4, y - ImGui::GetTextLineHeight() * 0.5f), IM_COL32(170, 170, 170, 200), buf);
    }
}

/**
 * 绘制 K 线（蜡烛图）。
 *
 * 参数说明：
 * - data: K 线数据序列（按时间升序）
 * - vs: 视图状态，包含横向缩放（scale_x）与滚动偏移（scroll_x）
 * - canvas_pos/canvas_size: 主图绘制区域左上角与尺寸（像素）
 * - begin/end: 需要绘制的数据索引区间 [begin, end)
 * - y_min/y_max: 映射到屏幕坐标的数值范围
 *
 * 说明：
 * - 上涨（收≥开）用绿色，下跌用红色。
 * - 细实体（|y_open−y_close| < 1 像素）退化为一条实线以提升可读性。
 * - 不负责裁剪，调用方需在外部 PushClipRect。
 */
static void draw_candles(const std::vector<Candle> &data, const ViewState &vs, const ImVec2 &canvas_pos, const ImVec2 &canvas_size,
                         int begin, int end, float y_min, float y_max) {
    ImDrawList *dl = ImGui::GetWindowDrawList();
    const float W = canvas_size.x;
    const float H = canvas_size.y;
    const float candle_w = std::max(1.0f, vs.scale_x * 0.7f);

    for (int i = begin; i < end; ++i) {
        const Candle &c = data[i];
        float x = canvas_pos.x + (i - vs.scroll_x) * vs.scale_x;
        float x0 = x - candle_w * 0.5f;
        float x1 = x + candle_w * 0.5f;
        auto y_to_screen = [&](double y) {
            float ty = (float)((y - y_min) / (y_max - y_min));
            return canvas_pos.y + (1.0f - ty) * H;
        };
        float y_open = y_to_screen(c.open);
        float y_close = y_to_screen(c.close);
        float y_high = y_to_screen(c.high);
        float y_low = y_to_screen(c.low);
        bool up = c.close >= c.open;
        ImU32 col = up ? IM_COL32(82, 196, 26, 255) : IM_COL32(255, 77, 79, 255);
        // Wick
        dl->AddLine(ImVec2(x, y_high), ImVec2(x, y_low), col, 1.0f);
        // Body
        if (std::abs(y_open - y_close) < 1.0f)
            dl->AddLine(ImVec2(x0, y_open), ImVec2(x1, y_close), col, candle_w);
        else
            dl->AddRectFilled(ImVec2(x0, std::min(y_open, y_close)), ImVec2(x1, std::max(y_open, y_close)), col);
    }
}

/**
 * 绘制一条按数据点连线的折线（指标/均线等）。
 *
 * 参数说明：
 * - series: 待绘制的数值序列（NAN 值会打断连线）
 * - vs: 视图状态（scale_x/scroll_x）
 * - canvas_pos/canvas_size: 绘制区域（像素）
 * - begin/end: 绘制索引范围 [begin, end)
 * - y_min/y_max: 垂直方向映射范围
 * - col: 线条颜色（ImU32）
 * - thickness: 线宽（像素）
 *
 * 说明：
 * - series[i] 为 NAN 时不连线，可用于产生断点。
 * - 仅负责连线，不包含裁剪；调用方需先设置 PushClipRect。
 */
static void draw_line_series(const std::vector<double> &series, const ViewState &vs,
                             const ImVec2 &canvas_pos, const ImVec2 &canvas_size,
                             int begin, int end, float y_min, float y_max, ImU32 col,
                             float thickness = 1.5f) {
    ImDrawList *dl = ImGui::GetWindowDrawList();
    const float H = canvas_size.y;
    auto y_to_screen = [&](double y) {
        float ty = (float)((y - y_min) / (y_max - y_min));
        return canvas_pos.y + (1.0f - ty) * H;
    };
    ImVec2 prev;
    bool has_prev = false;
    for (int i = begin; i < end; ++i) {
        double v = series[i];
        if (std::isnan(v)) {
            has_prev = false;
            continue;
        }
        float x = canvas_pos.x + (i - vs.scroll_x) * vs.scale_x;
        float y = y_to_screen(v);
        ImVec2 cur(x, y);
    if (has_prev) dl->AddLine(prev, cur, col, thickness);
        prev = cur;
        has_prev = true;
    }
}

// Compute Bollinger Bands (mid = SMA(period), upper/lower = mid ± k * stddev)
static void compute_boll(const std::vector<double>& series, int period, double k,
                         std::vector<double>& mid, std::vector<double>& upper, std::vector<double>& lower)
{
    size_t n = series.size();
    double NaN = std::numeric_limits<double>::quiet_NaN();
    mid.assign(n, NaN);
    upper.assign(n, NaN);
    lower.assign(n, NaN);
    if (period <= 1 || (size_t)period > n) return;
    for (size_t i = period - 1; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < period; ++j) sum += series[i - j];
        double mean = sum / (double)period;
        double var = 0.0;
        for (int j = 0; j < period; ++j) {
            double d = series[i - j] - mean;
            var += d * d;
        }
        var /= (double)period; // population variance
        double sd = std::sqrt(var);
        mid[i] = mean;
        upper[i] = mean + k * sd;
        lower[i] = mean - k * sd;
    }
}

void jump_window(const std::vector<Candle> &candles, ViewState &vs, const ImVec2 &main_size) {
    // Jump window: jump to date/time or to oldest/newest
    ImGui::Begin("Jump");
    static char date_buf[16] = ""; // YYYY-MM-DD
    static char time_buf[8] = "";  // HH:MM
    // prefill with  latest candle if empty
    if (date_buf[0] == '\0' || time_buf[0] == '\0') {
        int seed_idx = (int)candles.size() - 1;
        seed_idx = std::clamp(seed_idx, 0, (int)candles.size() - 1);
        if (!candles.empty()) {
            char tmp[64];
            format_time_label(candles[seed_idx].time, tmp, sizeof(tmp), true);
            // split into date and time
            // expected format: YYYY-MM-DD HH:MM
            const char *sp = std::strchr(tmp, ' ');
            if (sp) {
                std::snprintf(date_buf, sizeof(date_buf), "%.*s", (int)(sp - tmp), tmp);
                std::snprintf(time_buf, sizeof(time_buf), "%s", sp + 1);
            }
        }
    }
    ImGui::InputText("Date (YYYY-MM-DD)", date_buf, sizeof(date_buf));
    // ImGui::SameLine();
    ImGui::InputText("Time (HH:MM)", time_buf, sizeof(time_buf));
    bool jump_clicked = ImGui::Button("Jump");
    ImGui::SameLine();
    bool first_clicked = ImGui::Button("<< Oldest");
    ImGui::SameLine();
    bool last_clicked = ImGui::Button("Newest >>");
    int jump_to_idx = -1;
    if (jump_clicked) {
        std::time_t tt{};
        if (parse_datetime_to_time_t(date_buf, time_buf, tt)) {
            jump_to_idx = find_nearest_index_by_time(candles, tt);
        }
    }
    if (first_clicked && !candles.empty()) jump_to_idx = 0;
    if (last_clicked && !candles.empty()) jump_to_idx = (int)candles.size() - 1;
    if (jump_to_idx >= 0) {
        // Center the selected index in the main plot width
        float visible_count = std::max(1.0f, main_size.x / std::max(1.0f, vs.scale_x));
        float center_offset = visible_count * 0.5f;
        float upper = std::max(0.0f, (float)candles.size() - visible_count);
        vs.scroll_x = std::clamp((float)jump_to_idx - center_offset, 0.0f, upper);
    }
    ImGui::End();
}

int main(int, char **) {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    #if defined(__APPLE__)
    #define GL_SILENCE_DEPRECATION
    const char *glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    // #elif defined(__linux__)
    #else
    const char *glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    #endif

    GLFWwindow *window = glfwCreateWindow(1280, 800, "ImGui K-Line Viewer", NULL, NULL);
    if (window == NULL)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void)io;
    // Docking is not available on master branch by default
    ImGui::StyleColorsDark();
    // ImGui::StyleColorsLight(); // 蜡烛图需要自行更改浅色

    // Load Chinese font for UI (try a few relative paths)
    {
        const auto font_size = 15.0f;
        const char *font_rel1 = "微软雅黑.ttf";
        const char *font_rel2 = "../微软雅黑.ttf";
        const char *font_rel3 = "../../微软雅黑.ttf";
        ImFontConfig cfg;
        cfg.OversampleH = 2;
        cfg.OversampleV = 1;
        cfg.PixelSnapH = true;
        const ImWchar *ranges = io.Fonts->GetGlyphRangesChineseSimplifiedCommon();
        ImFont *font_cn = io.Fonts->AddFontFromFileTTF(font_rel1, font_size, &cfg, ranges);
        if (!font_cn) font_cn = io.Fonts->AddFontFromFileTTF(font_rel2, font_size, &cfg, ranges);
        if (!font_cn) font_cn = io.Fonts->AddFontFromFileTTF(font_rel3, font_size, &cfg, ranges);
        if (font_cn) io.FontDefault = font_cn; // use as default if loaded
    }

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Data
    std::vector<Candle> candles;      // active timeframe data for rendering
    std::vector<Candle> base_1m;      // sqlite 1m base when DB source is used
    // Try to load from CSV (several relative paths tried), fallback to synthetic
    const char *csv_rel1 = "test_data/sh000001(上证指数).csv";
    const char *csv_rel2 = "../test_data/sh000001(上证指数).csv";
    const char *csv_rel3 = "../../test_data/sh000001(上证指数).csv";
    std::string active_csv;
    if (load_csv_sh_index(csv_rel1, candles))
        active_csv = csv_rel1;
    else if (load_csv_sh_index(csv_rel2, candles))
        active_csv = csv_rel2;
    else if (load_csv_sh_index(csv_rel3, candles))
        active_csv = csv_rel3;
    else {
        candles = gen_sample_data(800);
        active_csv = "<synthetic>";
    }
    std::vector<double> closes;
    std::vector<double> highs;
    std::vector<double> lows;
    closes.reserve(candles.size());
    highs.reserve(candles.size());
    lows.reserve(candles.size());
    for (auto &c : candles) { closes.push_back(c.close); highs.push_back(c.high); lows.push_back(c.low); }

    // Indicator params (customizable)
    int sma_period = 20;
    int ema_period = 50;
    int macd_fast = 12, macd_slow = 26, macd_signal = 9;
    int rsi_period = 14;
    int kdj_period = 9;
    int boll_period = 20; // Bollinger period
    float boll_k = 2.0f;  // Bollinger multiplier
    // SAR & TD9 params
    float sar_step = 0.02f;
    float sar_max = 0.2f;
    int td_lookback = 4; // classic TD setup compares with 4 bars earlier

    // Indicator series
    std::vector<double> sma_v = ind::sma(closes, sma_period);
    std::vector<double> ema_v = ind::ema(closes, ema_period);
    std::vector<double> macd_line, signal_line, hist;
    ind::macd(closes, macd_fast, macd_slow, macd_signal, &macd_line, &signal_line, &hist);
    std::vector<double> rsi_v;
    // KDJ series
    std::vector<double> k_vals, d_vals, j_vals;
    ind::rsi(closes, rsi_period, rsi_v);
    ind::kdj(closes, highs, lows, kdj_period, k_vals, d_vals, j_vals);

    // Bollinger Bands series
    std::vector<double> boll_mid, boll_up, boll_dn;
    compute_boll(closes, boll_period, (double)boll_k, boll_mid, boll_up, boll_dn);
    // SAR & TD9 series
    std::vector<double> sar_v;
    std::vector<int> sar_trend; // +1 up, -1 down
    ind::parabolic_sar(highs, lows, closes, sar_step, sar_max, sar_v, &sar_trend);
    std::vector<int> td_buy, td_sell;
    ind::td_setup(closes, td_lookback, td_buy, td_sell);

    ViewState vs;
    ChartOptions opt;
    bool crosshair_visible = false;
    int cross_idx = -1;
    double cross_price = 0.0;
    // Trades state
    std::vector<TradeAnnot> trades;
    bool show_trades = true;

    // Data source control state
    enum class SourceType { CSV, SQLITE_1M };
    SourceType source = SourceType::CSV;
    // Timeframe list (seconds)
    struct Tf { const char* name; int secs; };
    static const Tf TF_LIST[] = {
        {"1m", 60}, {"5m", 300}, {"15m", 900}, {"30m", 1800}, {"1h", 3600}, {"4h", 14400}, {"1d", 86400}
    };
    int tf_index = 0; // default 1m for DB; for CSV it will be informational only
    // DB inputs
    DbLoadOptions dbopt;
    dbopt.db_path = "test_data/BTCUSDC.db";
    dbopt.symbol = "BTCUSDC";
    dbopt.start_ms = 0; // load all by default
    dbopt.end_ms = (int64_t)9e18;

    auto recompute_series = [&]() {
        closes.clear(); highs.clear(); lows.clear();
        closes.reserve(candles.size()); highs.reserve(candles.size()); lows.reserve(candles.size());
        for (auto &c : candles) { closes.push_back(c.close); highs.push_back(c.high); lows.push_back(c.low); }
        sma_v = ind::sma(closes, sma_period);
        ema_v = ind::ema(closes, ema_period);
        ind::macd(closes, macd_fast, macd_slow, macd_signal, &macd_line, &signal_line, &hist);
        ind::rsi(closes, rsi_period, rsi_v);
        compute_boll(closes, boll_period, (double)boll_k, boll_mid, boll_up, boll_dn);
    ind::parabolic_sar(highs, lows, closes, sar_step, sar_max, sar_v, &sar_trend);
    ind::td_setup(closes, td_lookback, td_buy, td_sell);
        // Reset scroll to show last if dataset changed significantly
        float visible_count = std::max(1.0f, 1200.0f / std::max(1.0f, vs.scale_x));
        if ((int)visible_count < (int)candles.size())
            vs.scroll_x = std::max(0.0f, (float)candles.size() - visible_count);
        else
            vs.scroll_x = 0.0f;
    };

    while (!glfwWindowShouldClose(window)) {
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

        // Reserve a right-side margin for floating labels so they don't overlap the plots
        static float right_margin = 110.0f; // px, adjustable in Legend
        float plot_width = std::max(10.0f, canvas_size.x - right_margin);
        float margin_x0 = canvas_pos.x + plot_width;
        float margin_x1 = canvas_pos.x + canvas_size.x;

    // Layout: main, optional RSI, optional MACD stacked vertically, plus a dedicated bottom time axis area
        const float spacing = 6.0f;
    const float axis_h = 28.0f; // dedicated bottom time axis height
        const float vol_h = opt.show_volume ? 80.0f : 0.0f;
        const float rsi_h = opt.show_rsi ? 80.0f : 0.0f;
        const float kdj_h = opt.show_kdj ? 80.0f : 0.0f;
        const float macd_h = opt.show_macd ? 120.0f : 0.0f;
    int gaps = 0;
        if (opt.show_volume) gaps++;
        if (opt.show_rsi) gaps++;
            if (opt.show_kdj) gaps++;
        if (opt.show_macd) gaps++;
    // add one extra gap above the time axis strip
    float main_h = canvas_size.y - axis_h - vol_h - rsi_h - kdj_h - macd_h - spacing * (gaps + 1);
        if (main_h < 120.0f) main_h = 120.0f;
        ImVec2 main_pos = canvas_pos;
        ImVec2 main_size = ImVec2(plot_width, main_h);
        ImVec2 vol_pos = ImVec2(main_pos.x, main_pos.y + main_size.y + (opt.show_volume ? spacing : 0.0f));
        ImVec2 vol_size = ImVec2(plot_width, vol_h);
        ImVec2 rsi_pos = opt.show_volume ? ImVec2(vol_pos.x, vol_pos.y + vol_size.y + (opt.show_rsi ? spacing : 0.0f)) : ImVec2(main_pos.x, main_pos.y + main_size.y + (opt.show_rsi ? spacing : 0.0f));
    ImVec2 rsi_size = ImVec2(plot_width, rsi_h);
    ImVec2 kdj_pos = (opt.show_rsi ? ImVec2(rsi_pos.x, rsi_pos.y + rsi_size.y + (opt.show_kdj ? spacing : 0.0f)) : (opt.show_volume ? ImVec2(vol_pos.x, vol_pos.y + vol_size.y + (opt.show_kdj ? spacing : 0.0f)) : ImVec2(main_pos.x, main_pos.y + main_size.y + (opt.show_kdj ? spacing : 0.0f))));
    ImVec2 kdj_size = ImVec2(plot_width, kdj_h);
        ImVec2 macd_pos = (opt.show_rsi ? ImVec2(rsi_pos.x, rsi_pos.y + rsi_size.y + (opt.show_macd ? spacing : 0.0f)) : (opt.show_volume ? ImVec2(vol_pos.x, vol_pos.y + vol_size.y + (opt.show_macd ? spacing : 0.0f)) : ImVec2(main_pos.x, main_pos.y + main_size.y + (opt.show_macd ? spacing : 0.0f))));
        ImVec2 macd_size = ImVec2(plot_width, macd_h);
    // Bottom time axis strip position/size (fixed height)
    ImVec2 axis_pos = ImVec2(main_pos.x, canvas_pos.y + canvas_size.y - axis_h);
    ImVec2 axis_size = ImVec2(plot_width, axis_h);


        // Interaction: mouse wheel zoom/scroll
        ImGuiIO &io = ImGui::GetIO();
        if (ImGui::IsWindowHovered()) {
            float wheel = io.MouseWheel;
            float wheel_h = io.MouseWheelH; // horizontal wheel (or shift+wheel) macos touchpad
            if (wheel != 0.0f && abs(wheel) > abs(wheel_h)) {
                float mouse_x = io.MousePos.x - main_pos.x;
                float center_index = vs.scroll_x + mouse_x / vs.scale_x;
                float prev_scale = vs.scale_x;
                vs.scale_x = std::clamp(vs.scale_x * (1.0f + wheel * 0.1f), 1.5f, 30.0f);
                vs.scroll_x = center_index - mouse_x / vs.scale_x;
            }
            if (wheel_h != 0.0f && abs(wheel) < abs(wheel_h)) {
                // move
                double c_move_num=1.5; //TODO
                vs.scroll_x -= c_move_num * wheel_h / vs.scale_x;
                vs.scroll_x = std::clamp(vs.scroll_x, 0.0f, (float)candles.size());
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
                bool in_rsi = opt.show_rsi && mp.x >= rsi_pos.x && mp.x <= rsi_pos.x + rsi_size.x && mp.y >= rsi_pos.y && mp.y <= rsi_pos.y + rsi_size.y;
                bool in_macd = opt.show_macd && mp.x >= macd_pos.x && mp.x <= macd_pos.x + macd_size.x && mp.y >= macd_pos.y && mp.y <= macd_pos.y + macd_size.y;
                if (in_main || in_rsi || in_macd) crosshair_visible = !crosshair_visible;
            }
        }

        // Determine visible range
        int begin = std::max(0, (int)std::floor(vs.scroll_x));
        int count = (int)std::ceil(main_size.x / vs.scale_x) + 2; // based on plot width only
        int end = std::min((int)candles.size(), begin + count);

        // Compute y-range for visible candles
        double y_min = 1e9, y_max = -1e9;
        for (int i = begin; i < end; ++i) {
            y_min = std::min(y_min, candles[i].low);
            y_max = std::max(y_max, candles[i].high);
        }
        float pad = (float)((y_max - y_min) * 0.05 + 1e-6);
        y_min -= pad;
        y_max += pad;

        // Draw backgrounds and borders per panel
        ImDrawList *dl = ImGui::GetWindowDrawList();
        auto rect = [&](ImVec2 p, ImVec2 s) { dl->AddRectFilled(p, ImVec2(p.x+s.x, p.y+s.y), IM_COL32(20,20,20,255)); dl->AddRect(p, ImVec2(p.x+s.x, p.y+s.y), IM_COL32(100,100,100,120)); };
        rect(main_pos, main_size);
        if (opt.show_volume) rect(vol_pos, vol_size);
    if (opt.show_volume) rect(vol_pos, vol_size);
    if (opt.show_rsi) dl->AddRect(rsi_pos, ImVec2(rsi_pos.x + rsi_size.x, rsi_pos.y + rsi_size.y), IM_COL32(100, 100, 100, 150));
    if (opt.show_kdj) dl->AddRect(kdj_pos, ImVec2(kdj_pos.x + kdj_size.x, kdj_pos.y + kdj_size.y), IM_COL32(100, 100, 100, 150));
    if (opt.show_macd) dl->AddRect(macd_pos, ImVec2(macd_pos.x + macd_size.x, macd_pos.y + macd_size.y), IM_COL32(100, 100, 100, 150));
        // KDJ panel
        if (opt.show_kdj) {
            dl->PushClipRect(kdj_pos, ImVec2(kdj_pos.x + kdj_size.x, kdj_pos.y + kdj_size.y), true);
            auto kdj_y = [&](double v) { return kdj_pos.y + (float)((100.0 - v) / 100.0) * kdj_size.y; };
            // 20/80 bands
            dl->AddLine(ImVec2(kdj_pos.x, kdj_y(20)), ImVec2(kdj_pos.x + kdj_size.x, kdj_y(20)), IM_COL32(150, 150, 150, 160));
            dl->AddLine(ImVec2(kdj_pos.x, kdj_y(80)), ImVec2(kdj_pos.x + kdj_size.x, kdj_y(80)), IM_COL32(150, 150, 150, 160));
            auto draw_line_local = [&](const std::vector<double>& s, ImU32 col) {
                ImVec2 prev; bool hp=false;
                for (int i = begin; i < end; ++i) {
                    double v = s[i];
                    if (std::isnan(v)) { hp=false; continue; }
                    float x = kdj_pos.x + (i - vs.scroll_x) * vs.scale_x;
                    float y = kdj_y(v);
                    ImVec2 cur(x,y);
                    if (hp) dl->AddLine(prev, cur, col, 1.5f);
                    prev = cur; hp=true;
                }
            };
            draw_line_local(k_vals, IM_COL32(64, 158, 255, 255));   // K - blue
            draw_line_local(d_vals, IM_COL32(255, 193, 7, 255));    // D - yellow
            draw_line_local(j_vals, IM_COL32(255, 99, 71, 255));    // J - tomato red
            dl->PopClipRect();
        }


        // Margin separator
        dl->AddLine(ImVec2(margin_x0, canvas_pos.y), ImVec2(margin_x0, canvas_pos.y + canvas_size.y), IM_COL32(100, 100, 100, 120));
    // Axis background (separate strip at bottom)
    dl->AddRectFilled(axis_pos, ImVec2(axis_pos.x + axis_size.x, axis_pos.y + axis_size.y), IM_COL32(18, 18, 18, 255));
    dl->AddRect(axis_pos, ImVec2(axis_pos.x + axis_size.x, axis_pos.y + axis_size.y), IM_COL32(100, 100, 100, 150));

    // Main chart (clip to plot area)
        dl->PushClipRect(main_pos, ImVec2(main_pos.x + main_size.x, main_pos.y + main_size.y), true);
        draw_grid(main_pos, ImVec2(main_pos.x + main_size.x, main_pos.y + main_size.y), (float)y_min, (float)y_max);
        draw_candles(candles, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max);
        if (opt.show_sma20) draw_line_series(sma_v, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, IM_COL32(255, 193, 7, 255));
        if (opt.show_ema50) draw_line_series(ema_v, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, IM_COL32(24, 144, 255, 255));
    if (opt.show_close_line) draw_line_series(closes, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, IM_COL32(220, 220, 220, 220));
        // SAR dots
        if (opt.show_sar) {
            const float H = main_size.y;
            auto y_to = [&](double y){ float ty=(float)((y - y_min)/(y_max - y_min)); return main_pos.y + (1.0f - ty) * H; };
            for (int i = begin; i < end; ++i) {
                double v = (i < (int)sar_v.size() ? sar_v[i] : NAN);
                if (std::isnan(v)) continue;
                float x = main_pos.x + (i - vs.scroll_x) * vs.scale_x;
                float y = y_to(v);
                int tr = (i < (int)sar_trend.size() ? sar_trend[i] : 0);
                ImU32 col = (tr >= 0 ? IM_COL32(64, 158, 255, 230) : IM_COL32(255, 99, 71, 230));
                ImU32 border = IM_COL32(0,0,0,200);
                dl->AddCircleFilled(ImVec2(x, y), opt.sar_dot_size, col);
                dl->AddCircle(ImVec2(x, y), opt.sar_dot_size, border, 0, 1.0f);
            }
        }
        // HLC area (Low-Close red fill, Close-High green fill) + three lines
        if (opt.show_hlc_area) {
            ImDrawList *dlx = dl;
            const float H = main_size.y;
            auto y_to = [&](double y){ float ty=(float)((y - y_min)/(y_max - y_min)); return main_pos.y + (1.0f - ty) * H; };
            ImU32 col_line_h = IM_COL32(82, 196, 26, 255);   // green high
            ImU32 col_line_c = IM_COL32(24, 144, 255, 255);  // blue close
            ImU32 col_line_l = IM_COL32(255, 77, 79, 255);   // red low
            ImU32 col_fill_up = IM_COL32(82, 196, 26, 60);   // translucent green
            ImU32 col_fill_dn = IM_COL32(255, 77, 79, 60);   // translucent red
            // Filled bands by quads between consecutive candles
            for (int i = begin; i < end - 1; ++i) {
                float x0 = main_pos.x + (i - vs.scroll_x) * vs.scale_x;
                float x1 = main_pos.x + ((i+1) - vs.scroll_x) * vs.scale_x;
                float yC0 = y_to(closes[i]);
                float yC1 = y_to(closes[i+1]);
                float yH0 = y_to(highs[i]);
                float yH1 = y_to(highs[i+1]);
                float yL0 = y_to(lows[i]);
                float yL1 = y_to(lows[i+1]);
                // Close -> High (green)
                ImVec2 gh[4] = { ImVec2(x0, yC0), ImVec2(x0, yH0), ImVec2(x1, yH1), ImVec2(x1, yC1) };
                dlx->AddConvexPolyFilled(gh, 4, col_fill_up);
                // Low -> Close (red)
                ImVec2 rl[4] = { ImVec2(x0, yL0), ImVec2(x0, yC0), ImVec2(x1, yC1), ImVec2(x1, yL1) };
                dlx->AddConvexPolyFilled(rl, 4, col_fill_dn);
            }
            // Lines: High (green), Close (blue), Low (red)
            draw_line_series(highs,  vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, col_line_h, 1.2f);
            draw_line_series(closes, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, col_line_c, 1.5f);
            draw_line_series(lows,   vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, col_line_l, 1.2f);
        }
        if (opt.show_boll) {
            ImU32 col_mid  = ImGui::ColorConvertFloat4ToU32(opt.boll_mid_color);
            ImU32 col_band = ImGui::ColorConvertFloat4ToU32(opt.boll_band_color);
            draw_line_series(boll_up,  vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, col_band, opt.boll_thickness);
            draw_line_series(boll_mid, vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, col_mid,  opt.boll_thickness);
            draw_line_series(boll_dn,  vs, main_pos, main_size, begin, end, (float)y_min, (float)y_max, col_band, opt.boll_thickness);
        }
        // TD9 numbers (setup counts 1..9)
        if (opt.show_td9) {
            auto y_to = [&](double y){ float ty=(float)((y - y_min)/(y_max - y_min)); return main_pos.y + (1.0f - ty) * main_size.y; };
            for (int i = begin; i < end; ++i) {
                int bc = (i < (int)td_buy.size() ? td_buy[i] : 0);
                int sc = (i < (int)td_sell.size() ? td_sell[i] : 0);
                if (bc > 0 && bc <= 9) {
                    char buf[8]; snprintf(buf, sizeof(buf), "%d", bc);
                    ImVec2 sz = ImGui::CalcTextSize(buf);
                    float x = main_pos.x + (i - vs.scroll_x) * vs.scale_x - sz.x * 0.5f;
                    float y = y_to(lows[i]) + 2.0f;
                    dl->AddText(ImVec2(x, y), IM_COL32(82, 196, 26, 240), buf);
                }
                if (sc > 0 && sc <= 9) {
                    char buf[8]; snprintf(buf, sizeof(buf), "%d", sc);
                    ImVec2 sz = ImGui::CalcTextSize(buf);
                    float x = main_pos.x + (i - vs.scroll_x) * vs.scale_x - sz.x * 0.5f;
                    float y = y_to(highs[i]) - sz.y - 2.0f;
                    dl->AddText(ImVec2(x, y), IM_COL32(255, 77, 79, 240), buf);
                }
            }
        }
        // Draw trades (markers + optional connecting line)
        if (show_trades && !trades.empty()) {
            auto x_from_time = [&](double t){ return main_pos.x + (float)((find_nearest_index_by_time(candles, (std::time_t)t) - vs.scroll_x) * vs.scale_x); };
            auto y_to = [&](double y){ float ty=(float)((y - y_min)/(y_max - y_min)); return main_pos.y + (1.0f - ty) * main_size.y; };
            for (auto &tr : trades) {
                float bx = x_from_time((double)tr.buy_t);
                float by = y_to(tr.buy_p);
                float sx = x_from_time((double)tr.sell_t);
                float sy = y_to(tr.sell_p);
                float r = 6.0f;
                // Buy marker: upward green triangle
                ImVec2 b0(bx, by - r), b1(bx - r, by + r), b2(bx + r, by + r);
                dl->AddTriangleFilled(b0, b1, b2, IM_COL32(82,196,26,220));
                dl->AddTriangle(b0, b1, b2, IM_COL32(0,0,0,200), 1.0f);
                // Sell marker: downward red triangle
                ImVec2 s0(sx, sy + r), s1(sx - r, sy - r), s2(sx + r, sy - r);
                dl->AddTriangleFilled(s0, s1, s2, IM_COL32(255,77,79,220));
                dl->AddTriangle(s0, s1, s2, IM_COL32(0,0,0,200), 1.0f);
                // Optional connecting line
                if (tr.connect_line) {
                    dl->AddLine(ImVec2(bx, by), ImVec2(sx, sy), IM_COL32(200,200,200,180), 1.5f);
                }
            }
        }
        dl->PopClipRect();

        // Price scale labels in the right margin for the main chart
        {
            dl->PushClipRect(ImVec2(margin_x0, main_pos.y), ImVec2(margin_x1, main_pos.y + main_size.y), true);
            int h_lines = 6; // match grid lines count
            double rng = (y_max - y_min);
            int decimals = (rng < 1.0 ? 4 : (rng < 10.0 ? 3 : 2));
            for (int i = 0; i <= h_lines; ++i) {
                float t = (float)i / (float)h_lines;
                float y = main_pos.y + (1.0f - t) * main_size.y;
                double price = y_max + (y_min - y_max) * t;
                char buf[32];
                if (decimals == 4)
                    snprintf(buf, sizeof(buf), "%.4f", price);
                else if (decimals == 3)
                    snprintf(buf, sizeof(buf), "%.3f", price);
                else
                    snprintf(buf, sizeof(buf), "%.2f", price);
                // small tick on the margin boundary
                dl->AddLine(ImVec2(margin_x0 + 1.0f, y), ImVec2(margin_x0 + 6.0f, y), IM_COL32(150, 150, 150, 160));
                ImVec2 ts = ImGui::CalcTextSize(buf);
                float ty = y - ts.y * 0.5f;
                dl->AddText(ImVec2(margin_x0 + 8.0f, ty), IM_COL32(200, 200, 200, 220), buf);
            }
            dl->PopClipRect();
        }

        // Volume bars
        if (opt.show_volume) {
            dl->PushClipRect(vol_pos, ImVec2(vol_pos.x + vol_size.x, vol_pos.y + vol_size.y), true);
            double vmax = 0.0;
            for (int i = begin; i < end; ++i) vmax = std::max(vmax, candles[i].volume);
            auto vy = [&](double v) { float t=(float)(v / (vmax + 1e-9)); return vol_pos.y + (1.0f - t) * vol_size.y; };
            float base_y = vol_pos.y + vol_size.y - 1.0f;
            for (int i = begin; i < end; ++i) {
                const auto &c = candles[i];
                float x = vol_pos.x + (i - vs.scroll_x) * vs.scale_x;
                float w = std::max(1.0f, vs.scale_x * 0.6f);
                float x0 = x - w * 0.5f, x1 = x + w * 0.5f;
                float y1 = vy(c.volume);
                ImU32 col = (c.close >= c.open) ? IM_COL32(82, 196, 26, 180) : IM_COL32(255, 77, 79, 180);
                dl->AddRectFilled(ImVec2(x0, base_y), ImVec2(x1, y1), col);
            }
            dl->PopClipRect();
        }

        // MACD panel
        if (opt.show_macd) {
            dl->PushClipRect(macd_pos, ImVec2(macd_pos.x + macd_size.x, macd_pos.y + macd_size.y), true);
            // Scale MACD
            double mmin = 1e9, mmax = -1e9;
            for (int i = begin; i < end; ++i) {
                mmin = std::min({mmin, macd_line[i], signal_line[i], hist[i]});
                mmax = std::max({mmax, macd_line[i], signal_line[i], hist[i]});
            }
            mmin = std::min(mmin, 0.0);
            mmax = std::max(mmax, 0.0);
            auto macd_y = [&](double v) { float t=(float)((v - mmin) / (mmax - mmin + 1e-9)); return macd_pos.y + (1.0f - t) * macd_size.y; };
            // Zero line
            float y0 = macd_y(0.0);
            dl->AddLine(ImVec2(macd_pos.x, y0), ImVec2(macd_pos.x + macd_size.x, y0), IM_COL32(120, 120, 120, 180));
            // Histogram
            for (int i = begin; i < end; ++i) {
                if (std::isnan(hist[i])) continue;
                float x = macd_pos.x + (i - vs.scroll_x) * vs.scale_x;
                float x0 = x - std::max(1.0f, vs.scale_x * 0.6f) / 2.0f;
                float x1 = x + std::max(1.0f, vs.scale_x * 0.6f) / 2.0f;
                float yv = macd_y(hist[i]);
                ImU32 col = hist[i] >= 0 ? IM_COL32(0, 200, 0, 180) : IM_COL32(220, 0, 0, 180);
                dl->AddRectFilled(ImVec2(x0, y0), ImVec2(x1, yv), col);
            }
            // Lines
            auto draw_line_local = [&](const std::vector<double> &s, ImU32 col) {
                ImVec2 p0 = macd_pos, sz = macd_size;
                const float H = sz.y;
                auto y_to = [&](double v) { float t=(float)((v - mmin) / (mmax - mmin + 1e-9)); return p0.y + (1.0f - t) * H; };
                ImVec2 prev;
                bool hp = false;
                for (int i = begin; i < end; ++i) {
                    double vv = s[i];
                    if (std::isnan(vv)) {
                        hp = false;
                        continue;
                    }
                    float x = p0.x + (i - vs.scroll_x) * vs.scale_x;
                    float y = y_to(vv);
                    ImVec2 cur(x, y);
                    if (hp) dl->AddLine(prev, cur, col, 1.5f);
                    prev = cur;
                    hp = true;
                }
            };
            draw_line_local(macd_line, IM_COL32(255, 255, 255, 220));
            draw_line_local(signal_line, IM_COL32(255, 215, 0, 220));
            dl->PopClipRect();
        }

        // RSI panel
        if (opt.show_rsi) {
            dl->PushClipRect(rsi_pos, ImVec2(rsi_pos.x + rsi_size.x, rsi_pos.y + rsi_size.y), true);
            auto rsi_y = [&](double v) { return rsi_pos.y + (float)((100.0 - v) / 100.0) * rsi_size.y; };
            // 30/70 bands
            dl->AddLine(ImVec2(rsi_pos.x, rsi_y(30)), ImVec2(rsi_pos.x + rsi_size.x, rsi_y(30)), IM_COL32(150, 150, 150, 180));
            dl->AddLine(ImVec2(rsi_pos.x, rsi_y(70)), ImVec2(rsi_pos.x + rsi_size.x, rsi_y(70)), IM_COL32(150, 150, 150, 180));
            // RSI line
            ImVec2 prev;
            bool has_prev = false;
            for (int i = begin; i < end; ++i) {
                double v = rsi_v[i];
                if (std::isnan(v)) {
                    has_prev = false;
                    continue;
                }
                float x = rsi_pos.x + (i - vs.scroll_x) * vs.scale_x;
                float y = rsi_y(v);
                ImVec2 cur(x, y);
                if (has_prev) dl->AddLine(prev, cur, IM_COL32(64, 158, 255, 255), 1.5f);
                prev = cur;
                has_prev = true;
            }
            dl->PopClipRect();
        }

        // Crosshair on panels and data readout
        auto format_opt = [](double v, char *buf, size_t n) { if (std::isnan(v)) { snprintf(buf,n,"-"); } else { snprintf(buf,n,"%.4f", v); } };
        if (crosshair_visible) {
            ImVec2 mp = io.MousePos;
            // Snap to nearest candle index based on mouse x
            float mouse_x_rel = std::clamp(mp.x - main_pos.x, 0.0f, main_size.x);
            int idx = (int)std::round(vs.scroll_x + mouse_x_rel / vs.scale_x);
            idx = std::clamp(idx, 0, (int)candles.size() - 1);
            cross_idx = idx;
            // Compute price at mouse y
            float mouse_y_rel = std::clamp(mp.y - main_pos.y, 0.0f, main_size.y);
            cross_price = (double)(y_max - (mouse_y_rel / main_size.y) * (y_max - y_min));

            // Crosshair lines
            float cx = main_pos.x + (cross_idx - vs.scroll_x) * vs.scale_x;
            float cy_main = main_pos.y + (float)((y_max - cross_price) / (y_max - y_min)) * main_size.y;
            dl->AddLine(ImVec2(cx, main_pos.y), ImVec2(cx, main_pos.y + main_size.y), IM_COL32(200, 200, 200, 120));
            dl->AddLine(ImVec2(main_pos.x, cy_main), ImVec2(main_pos.x + main_size.x, cy_main), IM_COL32(200, 200, 200, 120));
            if (opt.show_volume) dl->AddLine(ImVec2(cx, vol_pos.y), ImVec2(cx, vol_pos.y + vol_size.y), IM_COL32(200, 200, 200, 60));
            if (opt.show_rsi) dl->AddLine(ImVec2(cx, rsi_pos.y), ImVec2(cx, rsi_pos.y + rsi_size.y), IM_COL32(200, 200, 200, 60));
            if (opt.show_macd) dl->AddLine(ImVec2(cx, macd_pos.y), ImVec2(cx, macd_pos.y + macd_size.y), IM_COL32(200, 200, 200, 60));

            // Data box (colored indicators matching plot colors)
            const Candle &c = candles[cross_idx];
            char sv[32], ev[32], mv[32], sg[32], hs[32], rs[32];
            format_opt(sma_v[cross_idx], sv, sizeof(sv));
            format_opt(ema_v[cross_idx], ev, sizeof(ev));
            format_opt(macd_line[cross_idx], mv, sizeof(mv));
            format_opt(signal_line[cross_idx], sg, sizeof(sg));
            format_opt(hist[cross_idx], hs, sizeof(hs));
            format_opt(rsi_v[cross_idx], rs, sizeof(rs));
            char dt[64];
            format_time_label(c.time, dt, sizeof(dt), true);
            double chg_pct = NAN;
            double chg_abs = NAN;
            if (cross_idx > 0) {
                double pc = candles[cross_idx - 1].close;
                if (pc != 0.0) {
                    chg_abs = c.close - pc;
                    chg_pct = chg_abs / pc * 100.0;
                }
            }

            // Colors
            auto C = [](ImU32 u) { return ImGui::ColorConvertU32ToFloat4(u); };
            ImVec4 col_up = C(IM_COL32(82, 196, 26, 255));
            ImVec4 col_dn = C(IM_COL32(255, 77, 79, 255));
            ImVec4 col_sma = C(IM_COL32(255, 193, 7, 255));
            ImVec4 col_ema = C(IM_COL32(24, 144, 255, 255));
            ImVec4 col_rsi = C(IM_COL32(64, 158, 255, 255));
            ImVec4 col_macd = C(IM_COL32(255, 255, 255, 220));
            ImVec4 col_sig = C(IM_COL32(255, 215, 0, 220));
            ImVec4 col_hist_pos = C(IM_COL32(0, 200, 0, 220));
            ImVec4 col_hist_neg = C(IM_COL32(220, 0, 0, 220));

            ImVec2 box_pos = ImVec2(main_pos.x + 8, main_pos.y + 8);
            ImGui::SetCursorScreenPos(box_pos);
            ImGui::BeginChild("DataBox", ImVec2(460, 0), ImGuiChildFlags_Border | ImGuiChildFlags_AutoResizeY);
            ImGui::Text("Date: %s", dt);
            ImGui::Text("Idx: %d", cross_idx);
            ImGui::Text("O: %.4f  H: %.4f  L: %.4f  C: %.4f  V: %.2f", c.open, c.high, c.low, c.close, c.volume);
            if (!std::isnan(chg_abs) && !std::isnan(chg_pct)) {
                ImGui::TextColored(chg_abs >= 0 ? col_up : col_dn, "Chg: %s%.4f (%s%.2f%%)",
                                   (chg_abs >= 0 ? "+" : ""), chg_abs, (chg_pct >= 0 ? "+" : ""), chg_pct);
            } else {
                ImGui::Text("Chg: -");
            }
            // SMA / EMA
            ImGui::TextColored(col_sma, "SMA(%d): %s", sma_period, sv);
            ImGui::SameLine();
            ImGui::TextColored(col_ema, "  EMA(%d): %s", ema_period, ev);
            // MACD triple
            ImGui::Text("MACD(%d,%d,%d):", macd_fast, macd_slow, macd_signal);
            ImGui::SameLine();
            ImGui::TextColored(col_macd, "MACD %s", mv);
            ImGui::SameLine();
            ImGui::TextColored(col_sig, "  Sig %s", sg);
            {
                ImVec4 hcol = (hist[cross_idx] >= 0.0 || std::isnan(hist[cross_idx])) ? col_hist_pos : col_hist_neg;
                ImGui::SameLine();
                ImGui::TextColored(hcol, "  Hist %s", hs);
            }
            // RSI
            ImGui::TextColored(col_rsi, "RSI(%d): %s", rsi_period, rs);
            // KDJ in data box
            if (opt.show_kdj) {
                char kb[32], dbuf[32], jb[32];
                format_opt(k_vals[cross_idx], kb, sizeof(kb));
                format_opt(d_vals[cross_idx], dbuf, sizeof(dbuf));
                format_opt(j_vals[cross_idx], jb, sizeof(jb));
                ImGui::Text("KDJ(%d):", kdj_period);
                ImGui::SameLine(); ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(IM_COL32(64, 158, 255, 255)), "K %s", kb);
                ImGui::SameLine(); ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(IM_COL32(255, 193, 7, 255)),  " D %s", dbuf);
                ImGui::SameLine(); ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(IM_COL32(255, 99, 71, 255)),  " J %s", jb);
            }
            // BOLL values
            if (opt.show_boll) {
                auto C = [](ImVec4 v){ return v; };
                char bu[32], bm[32], bl[32];
                format_opt(boll_up[cross_idx], bu, sizeof(bu));
                format_opt(boll_mid[cross_idx], bm, sizeof(bm));
                format_opt(boll_dn[cross_idx], bl, sizeof(bl));
                ImVec4 col_mid  = opt.boll_mid_color;
                ImVec4 col_band = opt.boll_band_color;
                ImGui::Text("BOLL(%d,%.2f):", boll_period, boll_k);
                ImGui::SameLine(); ImGui::TextColored(C(col_mid),  "Mid %s", bm);
                ImGui::SameLine(); ImGui::TextColored(C(col_band), " Up %s", bu);
                ImGui::SameLine(); ImGui::TextColored(C(col_band), " Low %s", bl);
            }
            // SAR in data box
            if (opt.show_sar && cross_idx >= 0 && cross_idx < (int)sar_v.size()) {
                char sb[32];
                format_opt(sar_v[cross_idx], sb, sizeof(sb));
                int tr = (cross_idx < (int)sar_trend.size() ? sar_trend[cross_idx] : 0);
                ImVec4 scol = ImGui::ColorConvertU32ToFloat4(tr >= 0 ? IM_COL32(64, 158, 255, 230) : IM_COL32(255, 99, 71, 230));
                ImGui::TextColored(scol, "SAR: %s", sb);
            }
            // TD9 counts in data box
            if (opt.show_td9 && cross_idx >= 0 && cross_idx < (int)td_buy.size()) {
                int bc = td_buy[cross_idx];
                int sc = td_sell[cross_idx];
                ImGui::Text("TD9 L%d: ", td_lookback);
                ImGui::SameLine();
                ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(IM_COL32(82,196,26,240)), "B %d", bc);
                ImGui::SameLine();
                ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(IM_COL32(255,77,79,240)), "  S %d", sc);
            }
            ImGui::EndChild();
        }

        // Dedicated bottom time axis (ticks and labels)
        float axis_y = axis_pos.y + axis_size.y;
        dl->AddLine(ImVec2(axis_pos.x, axis_y), ImVec2(axis_pos.x + axis_size.x, axis_y), IM_COL32(120, 120, 120, 160));
        int px_step = 80; // min pixel between ticks
        int step = std::max(1, (int)std::round(px_step / std::max(1.0f, vs.scale_x)));
        int first = std::max(0, ((int)std::floor(vs.scroll_x) / step) * step);
        for (int i = first; i < (int)candles.size(); i += step) {
            float x = axis_pos.x + (i - vs.scroll_x) * vs.scale_x;
            if (x < axis_pos.x || x > axis_pos.x + axis_size.x) continue;
            dl->AddLine(ImVec2(x, axis_y), ImVec2(x, axis_y - 6.0f), IM_COL32(150, 150, 150, 160));
            char label[32];
            format_time_label(candles[i].time, label, sizeof(label), false);
            ImVec2 sz = ImGui::CalcTextSize(label);
            dl->AddText(ImVec2(x - sz.x * 0.5f, axis_y - sz.y - 8.0f), IM_COL32(180, 180, 180, 220), label);
        }
        // Crosshair time bubble in axis strip
        if (crosshair_visible && cross_idx >= 0 && cross_idx < (int)candles.size()) {
            float cx = axis_pos.x + (cross_idx - vs.scroll_x) * vs.scale_x;
            char tbuf[64];
            format_time_label(candles[cross_idx].time, tbuf, sizeof(tbuf), true);
            ImVec2 ts = ImGui::CalcTextSize(tbuf);
            float padx = 6.0f, pady = 3.0f;
            float w = ts.x + padx * 2.0f;
            float h = ts.y + pady * 2.0f;
            float x0 = std::clamp(cx - w * 0.5f, axis_pos.x, axis_pos.x + axis_size.x - w);
            float y0 = axis_pos.y + (axis_size.y - h) * 0.5f;
            ImVec2 p1(x0, y0), p2(x0 + w, y0 + h);
            dl->AddRectFilled(p1, p2, IM_COL32(40, 40, 40, 230), 4.0f);
            dl->AddRect(p1, p2, IM_COL32(0, 0, 0, 200), 4.0f, 0, 1.0f);
            dl->AddText(ImVec2(x0 + padx, y0 + pady), IM_COL32(240, 240, 240, 255), tbuf);
            // small top tick aligned with bubble center
            dl->AddLine(ImVec2(cx, axis_pos.y), ImVec2(cx, axis_pos.y + 6.0f), IM_COL32(180, 180, 180, 160));
        }

        // Right-side floating labels in the reserved margin for each panel (doesn't overlap plots)
        {
            // Compute the last truly visible candle by x-bound, so clipped candles are excluded
            int last = 0;
            if (!candles.empty()) {
                float candle_w = std::max(1.0f, vs.scale_x * 0.7f);
                float right_px = std::max(0.0f, main_size.x - candle_w * 0.5f);
                int last_by_x = (int)std::floor(vs.scroll_x + right_px / std::max(1.0f, vs.scale_x));
                last = std::clamp(last_by_x, begin, std::max(begin, end - 1));
            }
            // Choose label index: either follow crosshair (if enabled) or use last visible
            int label_idx = last;
            if (opt.labels_follow_cursor && crosshair_visible) {
                label_idx = std::clamp(cross_idx, begin, std::max(begin, end - 1));
            }
            if (label_idx >= begin && label_idx < end) {
                struct Label {
                    float y;
                    ImU32 col;
                    bool filled;
                    std::string text;
                };
                auto draw_labels = [&](const std::vector<Label> &in, ImVec2 panel_pos, ImVec2 panel_size) {
                    if (in.empty()) return;
                    std::vector<Label> labels = in;
                    std::sort(labels.begin(), labels.end(), [](const Label &a, const Label &b) { return a.y < b.y; });
                    float label_h = ImGui::GetTextLineHeightWithSpacing() + 6.0f;
                    for (size_t i = 1; i < labels.size(); ++i) {
                        if (labels[i].y - labels[i - 1].y < label_h) labels[i].y = labels[i - 1].y + label_h;
                    }
                    for (auto &L : labels) {
                        if (L.y < panel_pos.y) L.y = panel_pos.y;
                        if (L.y > panel_pos.y + panel_size.y - label_h) L.y = panel_pos.y + panel_size.y - label_h;
                    }
                    float x0 = margin_x0 + 4.0f, x1 = margin_x1 - 4.0f;
                    for (auto &L : labels) {
                        ImVec2 ts = ImGui::CalcTextSize(L.text.c_str());
                        float w = ts.x + 12.0f, h = ts.y + 6.0f;
                        ImVec2 p1 = ImVec2(x1 - w, L.y);
                        ImVec2 p2 = ImVec2(x1, L.y + h);
                        if (L.filled) {
                            dl->AddRectFilled(p1, p2, L.col, 4.0f);
                            dl->AddRect(p1, p2, IM_COL32(0, 0, 0, 180), 4.0f, 0, 1.0f);
                            dl->AddText(ImVec2(p1.x + 6.0f, p1.y + 3.0f), IM_COL32(255, 255, 255, 255), L.text.c_str());
                        } else {
                            dl->AddRect(p1, p2, L.col, 4.0f, 0, 1.5f);
                            dl->AddText(ImVec2(p1.x + 6.0f, p1.y + 3.0f), L.col, L.text.c_str());
                        }
                    }
                };

                // Main panel labels (price filled, SMA/EMA hollow)
                auto y_to_main = [&](double y) { float ty = (float)((y - y_min) / (y_max - y_min)); return main_pos.y + (1.0f - ty) * main_size.y; };
                std::vector<Label> main_labels;
                const Candle &lc = candles[label_idx];
                bool up = lc.close >= lc.open;
                char tbuf[64];
                snprintf(tbuf, sizeof(tbuf), "%.4f", lc.close);
                main_labels.push_back({y_to_main(lc.close), up ? IM_COL32(82, 196, 26, 255) : IM_COL32(255, 77, 79, 255), true, std::string(tbuf)});
                auto push_line = [&](const std::vector<double> &s, bool enabled, ImU32 col, const char *name) { if (!enabled) return; double v=s[label_idx]; if (std::isnan(v)) return; char b[64]; snprintf(b, sizeof(b), "%s %.4f", name, v); main_labels.push_back({ y_to_main(v), col, false, std::string(b) }); };
                push_line(sma_v, opt.show_sma20, IM_COL32(255, 193, 7, 255), "SMA");
                push_line(ema_v, opt.show_ema50, IM_COL32(24, 144, 255, 255), "EMA");
                // BOLL labels (hollow)
                if (opt.show_boll) {
                    ImU32 col_mid  = ImGui::ColorConvertFloat4ToU32(opt.boll_mid_color);
                    ImU32 col_band = ImGui::ColorConvertFloat4ToU32(opt.boll_band_color);
                    double vup = boll_up[label_idx];
                    double vm  = boll_mid[label_idx];
                    double vdn = boll_dn[label_idx];
                    if (!std::isnan(vup)) { char b[64]; snprintf(b, sizeof(b), "BOLL U %.4f", vup); main_labels.push_back({ y_to_main(vup), col_band, false, std::string(b) }); }
                    if (!std::isnan(vm))  { char b[64]; snprintf(b, sizeof(b), "BOLL M %.4f", vm);  main_labels.push_back({ y_to_main(vm),  col_mid,  false, std::string(b) }); }
                    if (!std::isnan(vdn)) { char b[64]; snprintf(b, sizeof(b), "BOLL L %.4f", vdn); main_labels.push_back({ y_to_main(vdn), col_band, false, std::string(b) }); }
                }
                draw_labels(main_labels, main_pos, main_size);

                // Volume labels (hollow)
                if (opt.show_volume) {
                    double vmax = 0.0;
                    for (int i = begin; i < end; ++i) vmax = std::max(vmax, candles[i].volume);
                    auto vy = [&](double v) { float t=(float)(v / (vmax + 1e-9)); return vol_pos.y + (1.0f - t) * vol_size.y; };
                    std::vector<Label> vol_labels;
                    char vb[64];
                    snprintf(vb, sizeof(vb), "VOL %.0f", candles[label_idx].volume);
                    vol_labels.push_back({vy(candles[label_idx].volume), IM_COL32(180, 180, 180, 220), false, std::string(vb)});
                    draw_labels(vol_labels, vol_pos, vol_size);
                }

                // RSI label (hollow, same color as RSI line)
                if (opt.show_rsi) {
                    auto rsi_y = [&](double v) { return rsi_pos.y + (float)((100.0 - v) / 100.0) * rsi_size.y; };
                    double rv = rsi_v[label_idx];
                    if (!std::isnan(rv)) {
                        std::vector<Label> rsi_labels;
                        char rb[64];
                        snprintf(rb, sizeof(rb), "RSI %.2f", rv);
                        rsi_labels.push_back({rsi_y(rv), IM_COL32(64, 158, 255, 255), false, std::string(rb)});
                        draw_labels(rsi_labels, rsi_pos, rsi_size);
                    }
                }

                // KDJ labels (hollow)
                if (opt.show_kdj) {
                    auto kdj_y = [&](double v) { return kdj_pos.y + (float)((100.0 - v) / 100.0) * kdj_size.y; };
                    std::vector<Label> kdj_labels;
                    auto pushk = [&](double v, ImU32 col, const char* nm){ if (!std::isnan(v)) { char b[64]; snprintf(b, sizeof(b), "%s %.2f", nm, v); kdj_labels.push_back({kdj_y(v), col, false, std::string(b)}); } };
                    pushk(k_vals[label_idx], IM_COL32(64, 158, 255, 255), "K");
                    pushk(d_vals[label_idx], IM_COL32(255, 193, 7, 255),  "D");
                    pushk(j_vals[label_idx], IM_COL32(255, 99, 71, 255),  "J");
                    draw_labels(kdj_labels, kdj_pos, kdj_size);
                }

                // MACD labels (hollow). Colors: line white, signal gold, hist green/red by sign
                if (opt.show_macd) {
                    double mmin = 1e9, mmax = -1e9;
                    for (int i = begin; i < end; ++i) {
                        mmin = std::min({mmin, macd_line[i], signal_line[i], hist[i]});
                        mmax = std::max({mmax, macd_line[i], signal_line[i], hist[i]});
                    }
                    mmin = std::min(mmin, 0.0);
                    mmax = std::max(mmax, 0.0);
                    auto macd_y = [&](double v) { float t=(float)((v - mmin) / (mmax - mmin + 1e-9)); return macd_pos.y + (1.0f - t) * macd_size.y; };
                    std::vector<Label> macd_labels;
                    double mv = macd_line[label_idx];
                    if (!std::isnan(mv)) {
                        char b[64];
                        snprintf(b, sizeof(b), "MACD %.4f", mv);
                        macd_labels.push_back({macd_y(mv), IM_COL32(255, 255, 255, 220), false, std::string(b)});
                    }
                    double sg = signal_line[label_idx];
                    if (!std::isnan(sg)) {
                        char b[64];
                        snprintf(b, sizeof(b), "Sig %.4f", sg);
                        macd_labels.push_back({macd_y(sg), IM_COL32(255, 215, 0, 220), false, std::string(b)});
                    }
                    double hv = hist[label_idx];
                    if (!std::isnan(hv)) {
                        ImU32 hc = hv >= 0 ? IM_COL32(0, 200, 0, 220) : IM_COL32(220, 0, 0, 220);
                        char b[64];
                        snprintf(b, sizeof(b), "Hist %.4f", hv);
                        macd_labels.push_back({macd_y(hv), hc, false, std::string(b)});
                    }
                    draw_labels(macd_labels, macd_pos, macd_size);
                }
            }
        }

        ImGui::End();

        jump_window(candles, vs, main_size);
        // Trades window: add trade by date/time and price
        ImGui::Begin("Trades");
        ImGui::Checkbox("Show trades", &show_trades);
        static char bdate[16] = "2025-01-02";
        static char btime[8]  = "09:30";
        static float bprice   = 3200.0f;
        static char sdate[16] = "2025-03-04";
        static char stime[8]  = "09:30";
        static float sprice   = 3300.0f;
        static bool  connect  = true;
        ImGui::InputText("Buy Date", bdate, sizeof(bdate));
        ImGui::InputText("Buy Time", btime, sizeof(btime));
        ImGui::InputFloat("Buy Price", &bprice, 0, 0, "%.4f");
        ImGui::Separator();
        ImGui::InputText("Sell Date", sdate, sizeof(sdate));
        ImGui::InputText("Sell Time", stime, sizeof(stime));
        ImGui::InputFloat("Sell Price", &sprice, 0, 0, "%.4f");
        ImGui::Checkbox("Connect line", &connect);
        if (ImGui::Button("Add Trade")) {
            std::time_t bt{}, st{};
            if (parse_datetime_to_time_t(bdate, btime, bt) && parse_datetime_to_time_t(sdate, stime, st)) {
                trades.push_back(TradeAnnot{bt, (double)bprice, st, (double)sprice, connect});
            }
        }
        // List trades
        for (int i = 0; i < (int)trades.size(); ++i) {
            ImGui::PushID(i);
            ImGui::Text("%d.", i+1);
            ImGui::SameLine(); ImGui::Text("B: %.4f S: %.4f", trades[i].buy_p, trades[i].sell_p);
            ImGui::SameLine();
            if (ImGui::Button("Remove")) { trades.erase(trades.begin()+i); ImGui::PopID(); break; }
            ImGui::PopID();
        }
        ImGui::End();

        // Separate Legend window
        ImGui::Begin("Legend");
        ImGui::Text("K-Line Viewer");
        ImGui::Text("Candles: %zu", candles.size());
        ImGui::Text("Data: %s", active_csv.c_str());
        ImGui::SliderFloat("Scale X", &vs.scale_x, 1.5f, 30.0f);
        ImGui::SliderFloat("Right margin", &right_margin, 60.0f, 240.0f, "%.0f px");
        ImGui::Text("Scroll: %.1f", vs.scroll_x);
        ImGui::Separator();
        // Data Source
        ImGui::Text("Data Source");
        int src = (source == SourceType::CSV ? 0 : 1);
        const char* src_names[] = {"CSV (SH Index)", "SQLite (1m)"};
        if (ImGui::Combo("Source", &src, src_names, IM_ARRAYSIZE(src_names))) {
            source = (src == 0 ? SourceType::CSV : SourceType::SQLITE_1M);
        }
        if (source == SourceType::SQLITE_1M) {
            static char db_path_buf[512] = "";
            static char symbol_buf[128] = "";
            if (db_path_buf[0] == '\0') std::snprintf(db_path_buf, sizeof(db_path_buf), "%s", dbopt.db_path.c_str());
            if (symbol_buf[0] == '\0') std::snprintf(symbol_buf, sizeof(symbol_buf), "%s", dbopt.symbol.c_str());
            ImGui::InputText("DB Path", db_path_buf, sizeof(db_path_buf));
            ImGui::InputText("Symbol", symbol_buf, sizeof(symbol_buf));
            ImGui::Text("Timeframe"); ImGui::SameLine();
            ImGui::Combo("##tf", &tf_index, [](void* data,int idx,const char** out_text){ auto list=(const Tf*)data; *out_text=list[idx].name; return true; }, (void*)TF_LIST, IM_ARRAYSIZE(TF_LIST));
            // Date/Time range inputs instead of milliseconds
            static char start_date[16] = ""; // YYYY-MM-DD
            static char start_time[8]  = ""; // HH:MM
            static char end_date[16]   = "";
            static char end_time[8]    = "";
            auto ms_to_dt = [](int64_t ms, char *d, size_t dn, char *t, size_t tn){
                if (ms <= 0) { if (d) d[0]='\0'; if (t) t[0]='\0'; return; }
                std::time_t tt = (std::time_t)(ms/1000);
                std::tm *lt = std::localtime(&tt);
                if (!lt) { if (d) d[0]='\0'; if (t) t[0]='\0'; return; }
                if (d) std::strftime(d, dn, "%Y-%m-%d", lt);
                if (t) std::strftime(t, tn, "%H:%M", lt);
            };
            // Prefill once
            if (start_date[0] == '\0' || start_time[0] == '\0') {
                if (!base_1m.empty()) {
                    int64_t ms = (int64_t)base_1m.front().time * 1000;
                    ms_to_dt(ms, start_date, sizeof(start_date), start_time, sizeof(start_time));
                } else if (dbopt.start_ms > 0) {
                    ms_to_dt(dbopt.start_ms, start_date, sizeof(start_date), start_time, sizeof(start_time));
                } else {
                    std::snprintf(start_date, sizeof(start_date), "2024-01-01");
                    std::snprintf(start_time, sizeof(start_time), "00:00");
                }
            }
            if (end_date[0] == '\0' || end_time[0] == '\0') {
                if (!base_1m.empty()) {
                    int64_t ms = (int64_t)base_1m.back().time * 1000;
                    ms_to_dt(ms, end_date, sizeof(end_date), end_time, sizeof(end_time));
                } else if (dbopt.end_ms < (int64_t)9e18) {
                    ms_to_dt(dbopt.end_ms, end_date, sizeof(end_date), end_time, sizeof(end_time));
                } else {
                    std::snprintf(end_date, sizeof(end_date), "2030-01-01");
                    std::snprintf(end_time, sizeof(end_time), "00:00");
                }
            }
            ImGui::InputText("Start Date", start_date, sizeof(start_date));
            ImGui::InputText("Start Time", start_time, sizeof(start_time));
            ImGui::InputText("End Date", end_date, sizeof(end_date));
            ImGui::InputText("End Time", end_time, sizeof(end_time));
            bool pressed = ImGui::Button("Load from DB");
            if (pressed) {
                dbopt.db_path = db_path_buf;
                dbopt.symbol = symbol_buf;
                // parse date/time -> ms
                std::time_t t0{}; std::time_t t1{};
                bool ok0 = parse_datetime_to_time_t(start_date, start_time, t0);
                bool ok1 = parse_datetime_to_time_t(end_date, end_time, t1);
                dbopt.start_ms = ok0 ? ((int64_t)t0 * 1000) : 0;
                dbopt.end_ms   = ok1 ? ((int64_t)t1 * 1000) : (int64_t)9e18;
                std::vector<Candle> tmp1m;
                if (load_sqlite_1m(dbopt, tmp1m)) {
                    base_1m.swap(tmp1m);
                    // aggregate
                    candles = aggregate_timeframe(base_1m, TF_LIST[tf_index].secs);
                    active_csv = dbopt.db_path + ":" + dbopt.symbol;
                    recompute_series();
                }
            }
            // Re-aggregate on tf change live (if base loaded)
            static int prev_tf = tf_index;
            if (tf_index != prev_tf) {
                prev_tf = tf_index;
                if (!base_1m.empty()) {
                    candles = aggregate_timeframe(base_1m, TF_LIST[tf_index].secs);
                    recompute_series();
                }
            }
        } else {
            // CSV source: timeframe is fixed by file; show combo disabled
            int dummy_tf = 0;
            ImGui::BeginDisabled();
            ImGui::Combo("##tf_csv", &dummy_tf, [](void* data,int idx,const char** out_text){ auto list=(const Tf*)data; *out_text=list[idx].name; return true; }, (void*)TF_LIST, IM_ARRAYSIZE(TF_LIST));
            ImGui::EndDisabled();
        }
        ImGui::Separator();
        ImGui::Text("Indicators");
        ImGui::Checkbox("SMA20", &opt.show_sma20);
        ImGui::Checkbox("EMA50", &opt.show_ema50);
        ImGui::Checkbox("MACD", &opt.show_macd);
        ImGui::Checkbox("RSI", &opt.show_rsi);
    ImGui::Checkbox("KDJ", &opt.show_kdj);
    ImGui::Checkbox("SAR 抛物线", &opt.show_sar);
    ImGui::Checkbox("神奇九转 (TD9)", &opt.show_td9);
        ImGui::Checkbox("Volume", &opt.show_volume);
        ImGui::Checkbox("Close line", &opt.show_close_line);
    ImGui::Checkbox("BOLL", &opt.show_boll);
    ImGui::Checkbox("HLC 区域", &opt.show_hlc_area);
        ImGui::Separator();
        ImGui::Text("Labels");
        ImGui::Checkbox("右侧浮标跟随光标索引", &opt.labels_follow_cursor);
        ImGui::Separator();
        ImGui::Text("Crosshair: %s (L-Click to toggle)", crosshair_visible ? "ON" : "OFF");
        ImGui::Separator();
        ImGui::Text("Parameters");
        bool dirty = false;
        dirty |= ImGui::SliderInt("SMA period", &sma_period, 2, 200);
        dirty |= ImGui::SliderInt("EMA period", &ema_period, 2, 200);
        dirty |= ImGui::SliderInt("MACD fast", &macd_fast, 2, 100);
        dirty |= ImGui::SliderInt("MACD slow", &macd_slow, 3, 200);
        dirty |= ImGui::SliderInt("MACD signal", &macd_signal, 2, 100);
    dirty |= ImGui::SliderInt("RSI period", &rsi_period, 2, 200);
    dirty |= ImGui::SliderInt("KDJ period", &kdj_period, 2, 200);
    dirty |= ImGui::SliderInt("BOLL period", &boll_period, 2, 300);
    dirty |= ImGui::SliderFloat("BOLL k", &boll_k, 0.5f, 4.0f, "%.2f");
    // SAR & TD9 params
    dirty |= ImGui::SliderFloat("SAR step", &sar_step, 0.001f, 0.2f, "%.3f");
    dirty |= ImGui::SliderFloat("SAR max", &sar_max, 0.01f, 1.0f, "%.2f");
    ImGui::SliderFloat("SAR dot size", &opt.sar_dot_size, 1.5f, 6.0f, "%.1f");
    dirty |= ImGui::SliderInt("TD9 lookback", &td_lookback, 1, 9);
    // BOLL style (no recompute required)
    ImGui::ColorEdit4("BOLL Mid Color", &opt.boll_mid_color.x, ImGuiColorEditFlags_NoInputs);
    ImGui::ColorEdit4("BOLL Band Color", &opt.boll_band_color.x, ImGuiColorEditFlags_NoInputs);
    ImGui::SliderFloat("BOLL thickness", &opt.boll_thickness, 0.5f, 4.0f, "%.2f");
        if (dirty) {
            sma_period = std::min(sma_period, (int)closes.size());
            ema_period = std::min(ema_period, (int)closes.size());
            macd_fast = std::max(2, macd_fast);
            macd_slow = std::max(macd_fast + 1, macd_slow);
            macd_signal = std::max(2, macd_signal);
            rsi_period = std::min(std::max(2, rsi_period), (int)closes.size());
            boll_period = std::min(std::max(2, boll_period), (int)closes.size());
            sar_step = std::max(0.0001f, std::min(1.0f, sar_step));
            sar_max = std::max(sar_step, sar_max);
            sma_v = ind::sma(closes, sma_period);
            ema_v = ind::ema(closes, ema_period);
            ind::macd(closes, macd_fast, macd_slow, macd_signal, &macd_line, &signal_line, &hist);
            ind::rsi(closes, rsi_period, rsi_v);
            compute_boll(closes, boll_period, (double)boll_k, boll_mid, boll_up, boll_dn);
            ind::kdj(closes, highs, lows, kdj_period, k_vals, d_vals, j_vals);
            ind::parabolic_sar(highs, lows, closes, sar_step, sar_max, sar_v, &sar_trend);
            ind::td_setup(closes, td_lookback, td_buy, td_sell);
        }
        ImGui::Separator();
        ImGui::Text("Load CSV");
        static char file_buf[512] = "";
        if (file_buf[0] == '\0') std::snprintf(file_buf, sizeof(file_buf), "%s", active_csv.c_str());
        ImGui::InputText("##csvpath", file_buf, sizeof(file_buf));
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
            std::vector<Candle> tmp;
            if (load_csv_sh_index(file_buf, tmp)) {
                source = SourceType::CSV;
                base_1m.clear();
                candles.swap(tmp);
                active_csv = file_buf;
                recompute_series();
            }
        }

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
