package com.hartwig.hmftools.patientdb;

import java.util.List;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Table;
import org.jooq.UpdatableRecord;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.testcontainers.mysql.MySQLContainer;

public abstract class DatabaseTestBase
{
    protected static final String TEST_SAMPLE_ID = "test_sample";

    protected static MySQLContainer container;
    protected static DatabaseAccess databaseAccess;

    @BeforeClass
    public static void createDatabase() throws Exception
    {
        container = new MySQLContainer("mysql:8")
                .withInitScript("generate_database.sql");

        container.start();

        databaseAccess = new DatabaseAccess(
                container.getUsername(),
                container.getPassword(),
                container.getJdbcUrl()
        );
    }

    @AfterClass
    public static void closeDatabase()
    {
        databaseAccess.close();
        container.stop();
    }

    public static <R extends UpdatableRecord<R>> List<R> fetchTable(Table<R> table, Class<R> recordClass)
    {
        return databaseAccess.context()
                .selectFrom(table)
                .fetchInto(recordClass);
    }
}
