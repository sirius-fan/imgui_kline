
# 示例代码
conn = sqlite3.connect(self.db_path)
cursor = conn.cursor()
        
        # 创建1分钟K线表
cursor.execute('''
            CREATE TABLE IF NOT EXISTS klines_1m (
                open_time INTEGER PRIMARY KEY,
                symbol TEXT NOT NULL,
                close_time INTEGER NOT NULL,
                open_price REAL NOT NULL,
                high_price REAL NOT NULL,
                low_price REAL NOT NULL,
                close_price REAL NOT NULL,
                volume REAL NOT NULL,
                quote_volume REAL NOT NULL,
                count INTEGER NOT NULL,
                taker_buy_volume REAL NOT NULL,
                taker_buy_quote_volume REAL NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # 创建1小时K线表
cursor.execute('''
            CREATE TABLE IF NOT EXISTS klines_1h (
                open_time INTEGER PRIMARY KEY,
                symbol TEXT NOT NULL,
                close_time INTEGER NOT NULL,
                open_price REAL NOT NULL,
                high_price REAL NOT NULL,
                low_price REAL NOT NULL,
                close_price REAL NOT NULL,
                volume REAL NOT NULL,
                quote_volume REAL NOT NULL,
                count INTEGER NOT NULL,
                taker_buy_volume REAL NOT NULL,
                taker_buy_quote_volume REAL NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
# 创建索引以提高查询性能
cursor.execute('CREATE INDEX IF NOT EXISTS idx_1m_symbol_time ON klines_1m(symbol, open_time)')
cursor.execute('CREATE INDEX IF NOT EXISTS idx_1h_symbol_time ON klines_1h(symbol, open_time)')
        
conn.commit()
conn.close()
logger.info(f"SQLite数据库初始化完成: {self.db_path}")