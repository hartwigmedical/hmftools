package com.hartwig.hmftools.patientdb.amber;

import static org.junit.Assert.assertEquals;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.patientdb.DatabaseTestBase;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadAmberDataTest extends DatabaseTestBase
{
    @Test
    public void shouldGenerateNewPatientIdForDifferentAmberSamples()
    {
        AmberSample amberSample = randomAmberSampleWithId("amber-sample-1");
        AmberSample otherAmberSample = randomAmberSampleWithId("amber-sample-2");

        LoadAmberData.processSample(amberSample, databaseAccess);
        LoadAmberData.processSample(otherAmberSample, databaseAccess);

        Set<Integer> loadedPatientIds = databaseAccess.readAmberPatients().stream().map(AmberPatient::patientId).collect(Collectors.toSet());

        assertEquals(2, loadedPatientIds.size());
    }

    @Test
    public void shouldReusePatientIdForSamplesBelongingToSamePatient()
    {
        AmberSample amberSample = randomAmberSampleWithId("amber-sample-1");
        AmberSample samePatientAmberSample = ImmutableAmberSample.builder().from(amberSample).sampleId("amber-sample-2").build();

        LoadAmberData.processSample(amberSample, databaseAccess);
        LoadAmberData.processSample(samePatientAmberSample, databaseAccess);

        Set<Integer> loadedPatientIds = databaseAccess.readAmberPatients().stream().map(AmberPatient::patientId).collect(Collectors.toSet());

        assertEquals(1, loadedPatientIds.size());
    }

    @Test
    public void shouldCorrectlyLoadAmberDataConcurrentlyAcrossSessions() throws SQLException
    {
        int nSamplesPerConnection = 100;
        int nConnections = 10;

        List<DatabaseAccess> connections = new ArrayList<>();
        for(int i = 0; i < nConnections; i++)
        {
            DatabaseAccess dbAccess = new DatabaseAccess(container.getUsername(), container.getPassword(), container.getJdbcUrl());
            connections.add(dbAccess);
        }

        try
        {
            // Asynchronously insert n samples per connection.
            // The insertion per connections still happens synchronously.
            List<List<AmberSample>> toInsertPerConnection = IntStream.range(0, nConnections)
                    .mapToObj(i -> IntStream.range(0, nSamplesPerConnection).mapToObj(j -> randomAmberSampleWithId(i + "-" + j)).toList())
                    .toList();

            Executor executor = Executors.newFixedThreadPool(nConnections);

            List<CompletableFuture<Void>> jobs = new ArrayList<>();
            for(int i = 0; i < toInsertPerConnection.size(); i++)
            {
                DatabaseAccess dbAccess = connections.get(i);
                List<AmberSample> amberSamples = toInsertPerConnection.get(i);

                jobs.add(CompletableFuture.runAsync(() -> {
                    for(AmberSample amberSample : amberSamples)
                    {
                        LoadAmberData.processSample(amberSample, dbAccess);
                    }
                }, executor));
            }

            jobs.forEach(CompletableFuture::join);

            // Assert there is a unique patientId for each patient in the db, and none got lost due to race conditions
            Set<Integer> patientIds = databaseAccess.readAmberPatients().stream().map(AmberPatient::patientId).collect(Collectors.toSet());
            assertEquals(nSamplesPerConnection * nConnections, patientIds.size());
        }
        finally
        {
            connections.forEach(DatabaseAccess::close);
        }
    }

    /**
     * Generate an AmberSample with random byte entries.
     * Although this method could generate two identical samples by chance, the odds of this happening are 1 in (2^8)^100
     */
    private static AmberSample randomAmberSampleWithId(String sampleId)
    {
        byte[] amberSites = new byte[100];
        new Random().nextBytes(amberSites);

        return ImmutableAmberSample.builder().sampleId(sampleId).entries(amberSites).build();
    }

}


