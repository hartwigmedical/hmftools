package com.hartwig.hmftools.patientdb.amber;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFFileReader;

public class LoadAmberData
{
    private static final int DEFAULT_MIN_DEPTH = 10;
    private static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    private static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    private static final String AMBER_SNP_VCF = "amber_snp_vcf";
    private static final String SNPCHECK_VCF = "snpcheck_vcf";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        String tumorSample = configBuilder.getValue(SAMPLE);
        String amberSnpPath = configBuilder.getValue(AMBER_SNP_VCF);
        String mappingLoci = configBuilder.getValue(SNPCHECK_VCF);

        LOGGER.info("Loading mapping loci from {}", mappingLoci);
        final ListMultimap<Chromosome, AmberSite> mappingSites = AmberSitesFile.sites(mappingLoci);

        final GenomePositionSelector<AmberSite> selector = GenomePositionSelectorFactory.create(mappingSites);

        try(final DatabaseAccess dbAccess = databaseAccess(configBuilder);

        final VCFFileReader fileReader = new VCFFileReader(new File(amberSnpPath), false))
        {
            LOGGER.info("Loading vcf snp data from {}", amberSnpPath);
            final List<BaseDepth> baseDepths = fileReader.iterator()
                    .stream()
                    .map(BaseDepth::fromVariantContext)
                    .filter(x -> selector.select(x).isPresent())
                    .collect(Collectors.toList());

            final AmberSampleFactory amberSampleFactory = new AmberSampleFactory(
                    DEFAULT_MIN_DEPTH, DEFAULT_MIN_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);

            final AmberSample sample = amberSampleFactory.fromBaseDepth(tumorSample, baseDepths);

            processSample(sample, dbAccess);
        }

        LOGGER.info("Amber data loading complete");
    }

    public static void processSample(final AmberSample sample, final DatabaseAccess dbAccess)
    {
        LOGGER.info("Comparing with existing samples");
        final List<AmberSample> allSamples = dbAccess.readAmberSamples();
        final List<AmberMapping> sampleMappings = new ArrayList<>();
        for(AmberSample other : allSamples)
        {
            if(!other.sampleId().equals(sample.sampleId()))
            {
                final AmberMapping mapping = AmberMappingFactory.create(sample, other);
                if(mapping.likelihood() > 0.8)
                {
                    sampleMappings.add(mapping);
                }
            }
        }

        LOGGER.info("Sample {} matched with {} other samples", sample.sampleId(), sampleMappings.size());
        final List<AmberPatient> existingPatients = dbAccess.readAmberPatients();
        final AmberPatientFactory amberPatientFactory = new AmberPatientFactory(existingPatients, sampleMappings);
        final AmberPatient patient = amberPatientFactory.createPatient(sample);

        LOGGER.info("Writing data");
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
    }
}
