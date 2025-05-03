package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_COUNT;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILE;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILTERS;
import static com.hartwig.hmftools.pave.PaveConfig.PON_ARTEFACTS_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
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

    public final List<AnnotationData> Annotators;

    public final RefGenomeInterface RefGenome;
    public final Reportability ReportableClassifier;

    public ReferenceData(final PaveConfig config, final ConfigBuilder configBuilder)
    {
        GeneDataCache = new GeneDataCache(
                configBuilder.getValue(ENSEMBL_DATA_DIR), config.RefGenVersion,
                configBuilder.getValue(DRIVER_GENE_PANEL));

        if(!GeneDataCache.loadCache(config.OnlyCanonical, false))
        {
            PV_LOGGER.error("gene data cache loading failed, exiting");
            System.exit(1);
        }

        ReportableClassifier = new Reportability(GeneDataCache.getDriverPanel());

        Annotators = Lists.newArrayList();

        // skip applying the PON and Gnomad annotations if already done by Sage
        VcfFileReader vcfFileReader = new VcfFileReader(config.VcfFile, true);
        boolean hasPonAnnotation = vcfFileReader.vcfHeader().hasInfoLine(PON_COUNT);
        boolean hasGnomadAnnotation = vcfFileReader.vcfHeader().hasInfoLine(GNOMAD_FREQ);

        Gnomad = new GnomadAnnotation(configBuilder, !hasGnomadAnnotation);
        Annotators.add(Gnomad);

        StandardPon = new PonAnnotation(!hasPonAnnotation ? configBuilder.getValue(PON_FILE) : null, true);
        StandardPon.loadFilters(configBuilder.getValue(PON_FILTERS));
        Annotators.add(StandardPon);

        ArtefactsPon = new PonAnnotation(configBuilder.getValue(PON_ARTEFACTS_FILE), false);
        Annotators.add(ArtefactsPon);

        VariantMappability = new Mappability(configBuilder);
        Annotators.add(VariantMappability);

        Clinvar = new ClinvarAnnotation(configBuilder);
        Annotators.add(Clinvar);

        BlacklistedVariants = new Blacklistings(configBuilder);
        Annotators.add(BlacklistedVariants);

        RefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
    }

    public boolean isValid()
    {
        boolean allValid = true;

        for(AnnotationData annotationData : Annotators)
        {
            if(annotationData.enabled() && !annotationData.hasValidData())
            {
                allValid = false;
                PV_LOGGER.error("invalid {} annotation data", annotationData.type());

            }
        }

        return allValid;
    }

    public void initialiseChromosomeData(final List<String> chromosomes, int threads)
    {
        final List<Callable> callableList = Lists.newArrayList();

        for(AnnotationData annotationData : Annotators)
        {
            if(annotationData.enabled())
            {
                annotationData.registerInitialChromosomes(chromosomes);

                if(annotationData instanceof Callable)
                    callableList.add((Callable)annotationData);
            }
        }

        TaskExecutor.executeTasks(callableList, threads);
    }

    public void onChromosomeComplete(final String chromosome)
    {
        for(AnnotationData annotationData : Annotators)
        {
            annotationData.onChromosomeComplete(chromosome);
        }
    }
}
