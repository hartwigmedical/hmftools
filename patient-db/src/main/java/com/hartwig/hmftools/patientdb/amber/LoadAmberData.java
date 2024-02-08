package com.hartwig.hmftools.patientdb.amber;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class LoadAmberData
{
    private static final int DEFAULT_MIN_DEPTH = 10;
    private static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    private static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;
    private static final double SAMPLE_MATCH_LIKELIHOOD = 0.8;

    private static final String AMBER_SNP_VCF = "amber_snp_vcf";
    private static final String SNPCHECK_VCF = "snpcheck_vcf";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        String tumorSample = configBuilder.getValue(SAMPLE);
        String amberSnpPath = configBuilder.getValue(AMBER_SNP_VCF);
        String snpCheckSitesVcf = configBuilder.getValue(SNPCHECK_VCF);

        final ListMultimap<Chromosome,AmberSite> sampleSnpCheckSites = AmberSitesFile.loadVcf(snpCheckSitesVcf);

        LOGGER.info("loaded {} SnpCheck sites from {}", sampleSnpCheckSites.size(), snpCheckSitesVcf);

        final AmberSampleFactory amberSampleFactory = new AmberSampleFactory(
                DEFAULT_MIN_DEPTH, DEFAULT_MIN_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);

        List<SiteEvidence> siteEvidenceList = Lists.newArrayListWithCapacity(sampleSnpCheckSites.size());

        DatabaseAccess dbAccess = databaseAccess(configBuilder);

        LOGGER.info("loading sample reference SnpCheck data from {}", amberSnpPath);

        VcfFileReader fileReader = new VcfFileReader(amberSnpPath);

        int refSiteCount = 0;
        for(VariantContext variant : fileReader.iterator())
        {
            SiteEvidence siteEvidence = SiteEvidence.fromVariantContext(variant);

            // check evidence sites vs required reference ones
            if(sampleSnpCheckSites.values().stream().anyMatch(x -> x.matches(
                    siteEvidence.Chromosome, siteEvidence.Position, siteEvidence.Ref, siteEvidence.Alt)))
            {
                ++refSiteCount;
                siteEvidenceList.add(siteEvidence);
            }
        }

        if(siteEvidenceList.size() != sampleSnpCheckSites.size())
        {
            LOGGER.error("sample({}) missing sites: required({}) matched({})",
                    tumorSample, sampleSnpCheckSites.size(), siteEvidenceList.size());
            System.exit(1);
        }

        LOGGER.debug("sample reference sites({}) filter to snpCheckSites({})", refSiteCount, siteEvidenceList.size());


        AmberSample sample = amberSampleFactory.createSampleData(tumorSample, siteEvidenceList);

        processSample(sample, dbAccess);

        LOGGER.info("Amber data loading complete");
    }

    public static void processSample(final AmberSample sample, final DatabaseAccess dbAccess)
    {
        LOGGER.info("comparing with existing samples");

        final List<AmberSample> allSamples = dbAccess.readAmberSamples();
        final List<AmberMapping> sampleMappings = new ArrayList<>();
        for(AmberSample other : allSamples)
        {
            if(!other.sampleId().equals(sample.sampleId()))
            {
                final AmberMapping mapping = AmberMappingFactory.create(sample, other);
                if(mapping.likelihood() > SAMPLE_MATCH_LIKELIHOOD)
                {
                    sampleMappings.add(mapping);
                }
            }
        }

        LOGGER.info("sample {} matched with {} other samples", sample.sampleId(), sampleMappings.size());

        final List<AmberPatient> existingPatients = dbAccess.readAmberPatients();
        final AmberPatientFactory amberPatientFactory = new AmberPatientFactory(existingPatients, sampleMappings);
        final AmberPatient patient = amberPatientFactory.createPatient(sample);

        LOGGER.info("writing sample data");
        dbAccess.writeAmberSample(sample);
        dbAccess.writeAmberMapping(sample.sampleId(), sampleMappings);
        dbAccess.writeAmberPatients(Collections.singletonList(patient));
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(AMBER_SNP_VCF, true, "Path to the amber snp vcf");
        configBuilder.addPath(SNPCHECK_VCF, true, "Path to the downsampled snp check vcf");
        addDatabaseCmdLineArgs(configBuilder, true);
        addLoggingOptions(configBuilder);
    }
}
