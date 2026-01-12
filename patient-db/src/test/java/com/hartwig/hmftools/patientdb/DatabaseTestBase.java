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
    protected static MySQLContainer CONTAINER;
    protected static DatabaseAccess DB_ACCESS;

    @BeforeClass
    public static void createDatabase() throws Exception
    {
        CONTAINER = new MySQLContainer("mysql:8")
                .withInitScript("generate_database.sql");

        CONTAINER.start();

        DB_ACCESS = new DatabaseAccess(
                CONTAINER.getUsername(),
                CONTAINER.getPassword(),
                CONTAINER.getJdbcUrl()
        );
    }

    @AfterClass
    public static void closeDatabase()
    {
        if(DB_ACCESS != null)
            DB_ACCESS.close();

        if(CONTAINER != null)
            CONTAINER.stop();
    }

    public static <R extends UpdatableRecord<R>> List<R> fetchTable(Table<R> table, Class<R> recordClass)
    {
        return DB_ACCESS.context()
                .selectFrom(table)
                .fetchInto(recordClass);
    }
}
