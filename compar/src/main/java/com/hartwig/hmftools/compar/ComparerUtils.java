package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CategoryType.DISRUPTION;
import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.chord.ChordComparer;
import com.hartwig.hmftools.compar.cider.Cdr3LocusSummaryComparer;
import com.hartwig.hmftools.compar.cider.CiderVdjComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.cuppa.CuppaComparer;
import com.hartwig.hmftools.compar.cuppa.CuppaImageComparer;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.isofox.IsofoxGeneDataComparer;
import com.hartwig.hmftools.compar.isofox.IsofoxSummaryComparer;
import com.hartwig.hmftools.compar.isofox.IsofoxTranscriptDataComparer;
import com.hartwig.hmftools.compar.isofox.NovelSpliceJunctionComparer;
import com.hartwig.hmftools.compar.isofox.RnaFusionComparer;
import com.hartwig.hmftools.compar.lilac.LilacAlleleComparer;
import com.hartwig.hmftools.compar.lilac.LilacQcComparer;
import com.hartwig.hmftools.compar.linx.DisruptionComparer;
import com.hartwig.hmftools.compar.linx.FusionComparer;
import com.hartwig.hmftools.compar.linx.GermlineSvComparer;
import com.hartwig.hmftools.compar.metrics.BamMetricsComparer;
import com.hartwig.hmftools.compar.metrics.FlagstatComparer;
import com.hartwig.hmftools.compar.mutation.GermlineVariantComparer;
import com.hartwig.hmftools.compar.mutation.SomaticVariantComparer;
import com.hartwig.hmftools.compar.peach.PeachComparer;
import com.hartwig.hmftools.compar.purple.CopyNumberComparer;
import com.hartwig.hmftools.compar.purple.GeneCopyNumberComparer;
import com.hartwig.hmftools.compar.purple.GermlineAmpDelComparer;
import com.hartwig.hmftools.compar.purple.PurityComparer;
import com.hartwig.hmftools.compar.sigs.SigsComparer;
import com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeComparer;
import com.hartwig.hmftools.compar.teal.TealComparer;
import com.hartwig.hmftools.compar.vchord.VChordComparer;
import com.hartwig.hmftools.compar.virus.VirusComparer;

public final class ComparerUtils
{
    public static List<ItemComparer> buildComparers(final ComparConfig config)
    {
        List<ItemComparer> comparers = Lists.newArrayList();
        MatchLevel matchLevel = config.MatchingLevel;

        // load in a predictable order irrespective of config
        for(CategoryType category : CategoryType.values())
        {
            if(!config.Categories.contains(category))
            {
                continue;
            }


            ItemComparer comparer = createComparer(category, config);

            if(matchLevel == REPORTABLE && !comparer.hasReportable())
            {
                continue;
            }

            comparers.add(comparer);
        }

        // link related or dependent comparers - could make this a virtual method too if becomes more common
        if(config.Categories.contains(FUSION))
        {
            FusionComparer fusionComparer = (FusionComparer)(comparers.stream()
                    .filter(x -> x.category() == FUSION).findFirst().orElse(null));

            DisruptionComparer disruptionComparer;
            if(config.Categories.contains(DISRUPTION))
            {
                disruptionComparer = (DisruptionComparer)(comparers.stream()
                        .filter(x -> x.category() == DISRUPTION).findFirst().orElse(null));
            }
            else
            {
                disruptionComparer = (DisruptionComparer)createComparer(DISRUPTION, config);
            }

            fusionComparer.setDisruptionComparer(disruptionComparer);
        }

        return comparers;
    }

    public static ItemComparer createComparer(final CategoryType category, final ComparConfig config)
    {
        switch(category)
        {
            case PURITY:
                return new PurityComparer(config);

            case DRIVER:
                return new DriverComparer(config);

            case COPY_NUMBER:
                return new CopyNumberComparer(config);

            case GENE_COPY_NUMBER:
                return new GeneCopyNumberComparer(config);

            case GERMLINE_AMP_DEL:
                return new GermlineAmpDelComparer(config);

            case FUSION:
                return new FusionComparer(config);

            case DISRUPTION:
                return new DisruptionComparer(config);

            case SOMATIC_VARIANT:
                return new SomaticVariantComparer(config);

            case GERMLINE_VARIANT:
                return new GermlineVariantComparer(config);

            case CUPPA:
                return new CuppaComparer(config);

            case CUPPA_IMAGE:
                return new CuppaImageComparer(config);

            case CHORD:
                return new ChordComparer(config);

            case LILAC_QC:
                return new LilacQcComparer(config);

            case LILAC_ALLELE:
                return new LilacAlleleComparer(config);

            case GERMLINE_SV:
                return new GermlineSvComparer(config);

            case PEACH:
                return new PeachComparer(config);

            case VIRUS:
                return new VirusComparer(config);

            case TUMOR_FLAGSTAT:
                return new FlagstatComparer(TUMOR_FLAGSTAT, config);

            case GERMLINE_FLAGSTAT:
                return new FlagstatComparer(GERMLINE_FLAGSTAT, config);

            case TUMOR_BAM_METRICS:
                return new BamMetricsComparer(TUMOR_BAM_METRICS, config);

            case GERMLINE_BAM_METRICS:
                return new BamMetricsComparer(GERMLINE_BAM_METRICS, config);

            case SNP_GENOTYPE:
                return new SnpGenotypeComparer(config);

            case CDR3_SEQUENCE:
                return new CiderVdjComparer(config);

            case CDR3_LOCUS_SUMMARY:
                return new Cdr3LocusSummaryComparer(config);

            case TELOMERE_LENGTH:
                return new TealComparer(config);

            case V_CHORD:
                return new VChordComparer(config);

            case SIGS:
                return new SigsComparer(config);

            case RNA_SUMMARY:
                return new IsofoxSummaryComparer(config);

            case RNA_GENE_DATA:
                return new IsofoxGeneDataComparer(config);

            case RNA_TRANSCRIPT_DATA:
                return new IsofoxTranscriptDataComparer(config);

            case NOVEL_SPLICE_JUNCTION:
                return new NovelSpliceJunctionComparer(config);

            case RNA_FUSION:
                return new RnaFusionComparer(config);

            default:
                return null;
        }
    }
}
