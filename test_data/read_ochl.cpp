#include <iostream>
#include <string>
#include <vector>
#include "sqlite3.h"

// 定义K线数据结构
struct KlineData {
    int64_t open_time;
    std::string symbol;
    int64_t close_time;
    double open_price;
    double high_price;
    double low_price;
    double close_price;
    double volume;
    double quote_volume;
    int count;
    double taker_buy_volume;
    double taker_buy_quote_volume;
    std::string created_at;
};

// 回调函数用于收集查询结果
static int klineCallback(void* data, int argc, char**argv, char** azColName) {
    if (data == nullptr) return 0;
    
    std::vector<KlineData>* klines = static_cast<std::vector<KlineData>*>(data);
    KlineData kline;
    
    // 解析查询结果到KlineData结构
    for (int i = 0; i < argc; i++) {
        std::string colName = azColName[i];
        std::string value = argv[i] ? argv[i] : "";
        
        if (colName == "open_time") kline.open_time = std::stoll(value);
        else if (colName == "symbol") kline.symbol = value;
        else if (colName == "close_time") kline.close_time = std::stoll(value);
        else if (colName == "open_price") kline.open_price = std::stod(value);
        else if (colName == "high_price") kline.high_price = std::stod(value);
        else if (colName == "low_price") kline.low_price = std::stod(value);
        else if (colName == "close_price") kline.close_price = std::stod(value);
        else if (colName == "volume") kline.volume = std::stod(value);
        else if (colName == "quote_volume") kline.quote_volume = std::stod(value);
        else if (colName == "count") kline.count = std::stoi(value);
        else if (colName == "taker_buy_volume") kline.taker_buy_volume = std::stod(value);
        else if (colName == "taker_buy_quote_volume") kline.taker_buy_quote_volume = std::stod(value);
        else if (colName == "created_at") kline.created_at = value;
    }
    
    klines->push_back(kline);
    return 0;
}

// 读取指定时间范围内的1分钟K线数据
bool read1mKlines(const std::string& dbPath, 
                 const std::string& symbol, 
                 int64_t startTime, 
                 int64_t endTime, 
                 std::vector<KlineData>& result) {
    sqlite3* db = nullptr;
    char* errMsg = nullptr;
    int rc;
    
    // 打开数据库
    rc = sqlite3_open(dbPath.c_str(), &db);
    if (rc != SQLITE_OK) {
        std::cerr << "无法打开数据库: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return false;
    }
    
    // 准备SQL查询语句
    std::string sql = "SELECT * FROM klines_1m "
                      "WHERE symbol = ? AND open_time >= ? AND open_time <= ? "
                      "ORDER BY open_time ASC;";
    
    sqlite3_stmt* stmt = nullptr;
    rc = sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr);
    
    if (rc != SQLITE_OK) {
        std::cerr << "SQL准备错误: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return false;
    }
    
    // 绑定参数
    sqlite3_bind_text(stmt, 1, symbol.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int64(stmt, 2, startTime);
    sqlite3_bind_int64(stmt, 3, endTime);
    
    // 执行查询并处理结果
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        KlineData kline;
        
        kline.open_time = sqlite3_column_int64(stmt, 0);
        kline.symbol = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        kline.close_time = sqlite3_column_int64(stmt, 2);
        kline.open_price = sqlite3_column_double(stmt, 3);
        kline.high_price = sqlite3_column_double(stmt, 4);
        kline.low_price = sqlite3_column_double(stmt, 5);
        kline.close_price = sqlite3_column_double(stmt, 6);
        kline.volume = sqlite3_column_double(stmt, 7);
        kline.quote_volume = sqlite3_column_double(stmt, 8);
        kline.count = sqlite3_column_int(stmt, 9);
        kline.taker_buy_volume = sqlite3_column_double(stmt, 10);
        kline.taker_buy_quote_volume = sqlite3_column_double(stmt, 11);
        kline.created_at = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 12));
        
        result.push_back(kline);
    }
    
    if (rc != SQLITE_DONE) {
        std::cerr << "查询执行错误: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_finalize(stmt);
        sqlite3_close(db);
        return false;
    }
    
    // 清理资源
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    return true;
}

// 打印K线数据
void printKlines(const std::vector<KlineData>& klines) {
    std::cout << "找到 " << klines.size() << " 条K线数据:" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "开盘时间\t交易对\t收盘时间\t开盘价\t最高价\t最低价\t收盘价\t成交量\t报价量\t成交笔数\t主动买入量\t主动买入报价量" << std::endl;
    
    for (const auto& kline : klines) {
        std::cout << kline.open_time << "\t"
                  << kline.symbol << "\t"
                  << kline.close_time << "\t"
                  << kline.open_price << "\t"
                  << kline.high_price << "\t"
                  << kline.low_price << "\t"
                  << kline.close_price << "\t"
                  << kline.volume << "\t"
                  << kline.quote_volume << "\t"
                  << kline.count << "\t"
                  << kline.taker_buy_volume << "\t"
                  << kline.taker_buy_quote_volume << std::endl;
    }
}

int main() {
    // 数据库路径
    std::string dbPath = "BTCUSDC.db";  // 替换为你的数据库路径
    
    // 查询参数
    std::string symbol = "BTCUSDC";          // 替换为你要查询的交易对
    int64_t startTime = 1727787600000;       // 起始时间戳(毫秒)
    int64_t endTime = 1755208800000;         // 结束时间戳(毫秒)
    
    // 存储查询结果
    std::vector<KlineData> klines;
    
    // 读取K线数据
    bool success = read1mKlines(dbPath, symbol, startTime, endTime, klines);
    
    if (success) {
        printKlines(klines);
    } else {
        std::cerr << "读取K线数据失败" << std::endl;
        return 1;
    }
    
    return 0;
}
// 输出示例
// 1755208620000   BTCUSDC 1755208679999   117912  117926  117900  117926  9.852   1.16165e+06     261   4.317    509024
// 1755208680000   BTCUSDC 1755208739999   117926  117933  117926  117928  7.55    890377  143     5.086 599794
// 1755208740000   BTCUSDC 1755208799999   117928  117942  117928  117942  5.619   662667  131     4.949 583655
// 1755208800000   BTCUSDC 1755208859999   117942  117942  117812  117812  17.722  2.08877e+06     393   1.604    189030