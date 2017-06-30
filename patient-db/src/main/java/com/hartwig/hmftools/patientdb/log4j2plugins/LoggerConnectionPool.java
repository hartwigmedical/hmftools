package com.hartwig.hmftools.patientdb.log4j2plugins;

import java.sql.Connection;
import java.sql.SQLException;

import javax.sql.DataSource;

import org.apache.commons.dbcp2.ConnectionFactory;
import org.apache.commons.dbcp2.DriverManagerConnectionFactory;
import org.apache.commons.dbcp2.PoolableConnection;
import org.apache.commons.dbcp2.PoolableConnectionFactory;
import org.apache.commons.dbcp2.PoolingDataSource;
import org.apache.commons.pool2.impl.GenericObjectPool;
import org.apache.commons.pool2.impl.GenericObjectPoolConfig;

public final class LoggerConnectionPool {
    private static final String dbUrl = "jdbc:mysql://localhost:3306/hmfpatients?useSSL=false&serverTimezone=CET";
    private static final String dbUser = "logger";
    private static final String dbPass = "Logger@1";

    private static final DataSource dataSource;

    private LoggerConnectionPool() {
    }

    static {
        final GenericObjectPoolConfig poolConfig = new GenericObjectPoolConfig();
        ConnectionFactory connectionFactory = new DriverManagerConnectionFactory(dbUrl, dbUser, dbPass);
        PoolableConnectionFactory poolableConnectionFactory = new PoolableConnectionFactory(connectionFactory, null);
        GenericObjectPool<PoolableConnection> pool = new GenericObjectPool<>(poolableConnectionFactory, poolConfig);
        poolableConnectionFactory.setPool(pool);
        dataSource = new PoolingDataSource<>(pool);
    }

    // MIVO: getDatabaseConnection is used by the JDBC log appender. see resources/log4j2.xml
    public static Connection getDatabaseConnection() throws SQLException {
        return dataSource.getConnection();
    }
}
