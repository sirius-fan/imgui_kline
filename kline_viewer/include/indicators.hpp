#pragma once
#include <vector>
#include <numeric>
#include <cmath>
#include <cstdint>
#include <algorithm>

struct Candle {
    uint64_t time;   // seconds or index
    double open;
    double high;
    double low;
    double close;
    double volume;
};

namespace ind {

inline std::vector<double> sma(const std::vector<double>& data, int period) {
    std::vector<double> out(data.size(), NAN);
    if (period <= 0 || (int)data.size() < period) return out;
    double sum = std::accumulate(data.begin(), data.begin()+period, 0.0);
    out[period-1] = sum / period;
    for (size_t i = period; i < data.size(); ++i) {
        sum += data[i] - data[i - period];
        out[i] = sum / period;
    }
    return out;
}

inline std::vector<double> ema(const std::vector<double>& data, int period) {
    std::vector<double> out(data.size(), NAN);
    if (period <= 0 || data.empty()) return out;
    double k = 2.0 / (period + 1.0);
    double ema_prev = data[0];
    out[0] = ema_prev;
    for (size_t i = 1; i < data.size(); ++i) {
        ema_prev = data[i] * k + ema_prev * (1.0 - k);
        out[i] = ema_prev;
    }
    return out;
}

inline void macd(const std::vector<double>& data, int fast=12, int slow=26, int signal=9,
                 std::vector<double>* macd_line=nullptr,
                 std::vector<double>* signal_line=nullptr,
                 std::vector<double>* histogram=nullptr) {
    auto ema_fast = ema(data, fast);
    auto ema_slow = ema(data, slow);
    std::vector<double> macd_v(data.size(), NAN);
    for (size_t i = 0; i < data.size(); ++i) macd_v[i] = ema_fast[i] - ema_slow[i];
    if (macd_line) *macd_line = macd_v;
    if (signal_line) *signal_line = ema(macd_v, signal);
    if (histogram) {
        auto sig = ema(macd_v, signal);
        histogram->resize(data.size(), NAN);
        for (size_t i = 0; i < data.size(); ++i) (*histogram)[i] = macd_v[i] - sig[i];
    }
}

inline void rsi(const std::vector<double>& data, int period, std::vector<double>& out) {
    out.assign(data.size(), NAN);
    if ((int)data.size() <= period) return;
    double gain = 0.0, loss = 0.0;
    for (int i = 1; i <= period; ++i) {
        double diff = data[i] - data[i-1];
        if (diff >= 0) gain += diff; else loss -= diff;
    }
    double avg_gain = gain / period;
    double avg_loss = loss / period;
    out[period] = 100.0 - 100.0 / (1.0 + (avg_loss == 0 ? 0 : (avg_gain/avg_loss)));
    for (size_t i = period+1; i < data.size(); ++i) {
        double diff = data[i] - data[i-1];
        double g = diff > 0 ? diff : 0.0;
        double l = diff < 0 ? -diff : 0.0;
        avg_gain = (avg_gain * (period - 1) + g) / period;
        avg_loss = (avg_loss * (period - 1) + l) / period;
        double rs = (avg_loss == 0) ? 0 : (avg_gain / avg_loss);
        out[i] = 100.0 - 100.0 / (1.0 + rs);
    }
}

// KDJ indicator (RSV period n, smoothing 1/3 for K and D)
inline void kdj(const std::vector<double>& close,
                const std::vector<double>& high,
                const std::vector<double>& low,
                int period,
                std::vector<double>& K,
                std::vector<double>& D,
                std::vector<double>& J) {
    size_t n = close.size();
    K.assign(n, NAN);
    D.assign(n, NAN);
    J.assign(n, NAN);
    if (period <= 0 || high.size() != n || low.size() != n || n == 0) return;
    double k_prev = 50.0, d_prev = 50.0; // common initialization
    for (size_t i = 0; i < n; ++i) {
        if (i + 1 < (size_t)period) {
            // not enough lookback yet, keep prev smoothed values
            K[i] = k_prev; D[i] = d_prev; J[i] = 3.0 * K[i] - 2.0 * D[i];
            continue;
        }
        size_t start = i + 1 - (size_t)period;
        double hh = high[start];
        double ll = low[start];
        for (size_t j = start + 1; j <= i; ++j) { hh = std::max(hh, high[j]); ll = std::min(ll, low[j]); }
        double rsv = 0.0;
        double denom = (hh - ll);
        if (denom == 0) rsv = 0.0;
        else rsv = (close[i] - ll) / denom * 100.0;
        double k = (2.0/3.0) * k_prev + (1.0/3.0) * rsv;
        double d = (2.0/3.0) * d_prev + (1.0/3.0) * k;
        double j = 3.0 * k - 2.0 * d;
        K[i] = k; D[i] = d; J[i] = j;
        k_prev = k; d_prev = d;
    }
}

// Parabolic SAR (Welles Wilder)
// step: acceleration factor step (e.g., 0.02), max_step: maximum AF (e.g., 0.2)
// out: SAR values; trend(optional): +1 for uptrend, -1 for downtrend
inline void parabolic_sar(const std::vector<double>& high,
                          const std::vector<double>& low,
                          const std::vector<double>& close,
                          double step,
                          double max_step,
                          std::vector<double>& out,
                          std::vector<int>* trend = nullptr) {
    size_t n = high.size();
    out.assign(n, NAN);
    if (trend) trend->assign(n, 0);
    if (n == 0 || low.size() != n || close.size() != n) return;
    // Initial direction by first two closes
    bool up = (n >= 2 ? (close[1] >= close[0]) : true);
    double af = step;
    double ep = up ? high[0] : low[0];
    double sar = up ? low[0] : high[0];
    for (size_t i = 0; i < n; ++i) {
        if (i == 0) {
            out[i] = sar;
            if (trend) (*trend)[i] = up ? 1 : -1;
            continue;
        }
        // Advance SAR
        sar = sar + af * (ep - sar);
        // Clamp to previous 1-2 bars extremes to avoid penetration
        if (up) {
            double ll1 = low[i-1];
            double ll2 = (i >= 2 ? low[i-2] : ll1);
            sar = std::min(sar, std::min(ll1, ll2));
        } else {
            double hh1 = high[i-1];
            double hh2 = (i >= 2 ? high[i-2] : hh1);
            sar = std::max(sar, std::max(hh1, hh2));
        }
        // Check reversal
        bool reversal = false;
        if (up) {
            if (low[i] < sar) { // reverse to downtrend
                up = false;
                sar = ep;               // on reversal, SAR set to prior EP
                af = step;              // reset AF
                ep = low[i];            // new EP for downtrend
                reversal = true;
            }
        } else {
            if (high[i] > sar) { // reverse to uptrend
                up = true;
                sar = ep;
                af = step;
                ep = high[i];
                reversal = true;
            }
        }
        if (!reversal) {
            // Update EP and AF within same trend
            if (up) {
                if (high[i] > ep) {
                    ep = high[i];
                    af = std::min(max_step, af + step);
                }
            } else {
                if (low[i] < ep) {
                    ep = low[i];
                    af = std::min(max_step, af + step);
                }
            }
        }
        out[i] = sar;
        if (trend) (*trend)[i] = up ? 1 : -1;
    }
}

// TD Sequential Setup ("神奇九转"基础版): counts vs close[i - lookback]
// Produces buy and sell setup counts (reset when condition fails). Counts grow 1..N (we'll draw up to 9).
inline void td_setup(const std::vector<double>& close,
                     int lookback,
                     std::vector<int>& buyCount,
                     std::vector<int>& sellCount) {
    size_t n = close.size();
    buyCount.assign(n, 0);
    sellCount.assign(n, 0);
    if (n == 0 || lookback <= 0) return;
    for (size_t i = 0; i < n; ++i) {
        if (i < (size_t)lookback) { buyCount[i] = 0; sellCount[i] = 0; continue; }
        // Buy setup: close < close[i-lookback]
        if (close[i] < close[i - lookback])
            buyCount[i] = (buyCount[i-1] > 0 ? buyCount[i-1] + 1 : 1);
        else
            buyCount[i] = 0;
        // Sell setup: close > close[i-lookback]
        if (close[i] > close[i - lookback])
            sellCount[i] = (sellCount[i-1] > 0 ? sellCount[i-1] + 1 : 1);
        else
            sellCount[i] = 0;
    }
}

} // namespace ind
