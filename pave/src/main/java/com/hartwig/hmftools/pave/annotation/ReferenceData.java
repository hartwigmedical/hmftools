package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.pave.PaveConfig.PON_ARTEFACTS_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PON_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PON_FILTERS;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.PaveConfig;

public class ReferenceData
{
    public final com.hartwig.hmftools.pave.GeneDataCache GeneDataCache;
    public final GnomadAnnotation Gnomad;
    public final PonAnnotation StandardPon;
    public final PonAnnotation ArtefactsPon;
    public final Mappability VariantMappability;
    public final ClinvarAnnotation Clinvar;
    public final Blacklistings BlacklistedVariants;
    public final RefGenomeInterface RefGenome;
    public final Reportability ReportableClassifier;

    public ReferenceData(final PaveConfig config, final ConfigBuilder configBuilder)
    {
        GeneDataCache = new GeneDataCache(
                configBuilder.getValue(ENSEMBL_DATA_DIR), config.RefGenVersion,
                configBuilder.getValue(DRIVER_GENE_PANEL_OPTION), true);

        if(!GeneDataCache.loadCache(config.OnlyCanonical, false))
        {
            PV_LOGGER.error("gene data cache loading failed, exiting");
            System.exit(1);
        }

        ReportableClassifier = new Reportability(GeneDataCache.getDriverPanel());

        Gnomad = new GnomadAnnotation(configBuilder);

        StandardPon = new PonAnnotation(configBuilder.getValue(PON_FILE), true);
        StandardPon.loadFilters(configBuilder.getValue(PON_FILTERS));

        ArtefactsPon = new PonAnnotation(configBuilder.getValue(PON_ARTEFACTS_FILE), false);

        VariantMappability = new Mappability(configBuilder);
        Clinvar = new ClinvarAnnotation(configBuilder);
        BlacklistedVariants = new Blacklistings(configBuilder);

        RefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
    }

    public boolean isValid()
    {
        if(StandardPon.enabled() && !StandardPon.hasValidData())
        {
            PV_LOGGER.error("invalid PON");
            return false;
        }

        if(ArtefactsPon.enabled() && !ArtefactsPon.hasValidData())
        {
            PV_LOGGER.error("invalid PON artefacts");
            return false;
        }

        if(!VariantMappability.hasValidData())
        {
            PV_LOGGER.error("invalid mappability data");
            return false;
        }

        if(!Clinvar.hasValidData())
        {
            PV_LOGGER.error("invalid Clinvar data");
            return false;
        }

        if(!BlacklistedVariants.hasValidData())
        {
            PV_LOGGER.error("invalid blacklistings data");
            return false;
        }

        if(!Gnomad.hasValidData())
        {
            PV_LOGGER.error("invalid Gnomad data");
            return false;
        }

        return true;
    }

    public void initialiseChromosomeData(final String chromosome)
    {
        if(Gnomad.hasData())
            Gnomad.getChromosomeCache(chromosome);

        if(VariantMappability.hasData())
            VariantMappability.getChromosomeCache(chromosome);

        if(StandardPon.enabled())
            StandardPon.getChromosomeCache(chromosome);
    }

    public void initialiseChromosomeData(final List<String> chromosomes, int threads)
    {
        final List<Callable> callableList = Lists.newArrayList();

        if(Gnomad.enabled())
        {
            Gnomad.registerInitialChromosomes(chromosomes);
            callableList.add(Gnomad);
        }

        if(VariantMappability.enabled())
        {
            VariantMappability.registerInitialChromosomes(chromosomes);
            callableList.add(VariantMappability);
        }

        if(StandardPon.enabled())
        {
            StandardPon.registerInitialChromosomes(chromosomes);
            callableList.add(StandardPon);
        }

        if(Clinvar.enabled())
        {
            callableList.add(Clinvar);
        }

        TaskExecutor.executeTasks(callableList, threads);
    }
}
