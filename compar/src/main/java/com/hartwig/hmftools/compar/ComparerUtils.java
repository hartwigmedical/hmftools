package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CategoryType.DISRUPTION;
import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.chord.ChordComparer;
import com.hartwig.hmftools.compar.cider.Cdr3LocusSummaryComparer;
import com.hartwig.hmftools.compar.cider.CiderVdjComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.cuppa.CuppaComparer;
import com.hartwig.hmftools.compar.cuppa.CuppaImageComparer;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.isofox.RnaGeneDataComparer;
import com.hartwig.hmftools.compar.isofox.RnaSummaryComparer;
import com.hartwig.hmftools.compar.isofox.RnaTranscriptDataComparer;
import com.hartwig.hmftools.compar.isofox.RnaNovelSpliceJunctionComparer;
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

import org.jetbrains.annotations.Nullable;

public final class ComparerUtils
{
    public static List<ItemComparer> buildComparers(final ComparConfig config, final FieldCheckCache fieldCheckCache)
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

            ItemComparer comparer = createComparer(category, config, fieldCheckCache);

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
                disruptionComparer = (DisruptionComparer)createComparer(DISRUPTION, config, fieldCheckCache);
            }

            fusionComparer.setDisruptionComparer(disruptionComparer);
        }

        return comparers;
    }

    public static ItemComparer createComparer(
            final CategoryType category, final ComparConfig config, @Nullable final FieldCheckCache fieldCheckCache)
    {
        Map<String,FieldCheck> categoryFieldOverrides =
                fieldCheckCache != null ? fieldCheckCache.getCategoryFieldCheckOverrides(category) : Collections.emptyMap();

        switch(category)
        {
            case PURITY:
                return new PurityComparer(config, categoryFieldOverrides);

            case DRIVER:
                return new DriverComparer(config, categoryFieldOverrides);

            case COPY_NUMBER:
                return new CopyNumberComparer(config, categoryFieldOverrides);

            case GENE_COPY_NUMBER:
                return new GeneCopyNumberComparer(config, categoryFieldOverrides);

            case GERMLINE_AMP_DEL:
                return new GermlineAmpDelComparer(config, categoryFieldOverrides);

            case FUSION:
                return new FusionComparer(config, categoryFieldOverrides);

            case DISRUPTION:
                return new DisruptionComparer(config, categoryFieldOverrides);

            case SOMATIC_VARIANT:
                return new SomaticVariantComparer(config, categoryFieldOverrides);

            case GERMLINE_VARIANT:
                return new GermlineVariantComparer(config, categoryFieldOverrides);

            case CUPPA:
                return new CuppaComparer(config, categoryFieldOverrides);

            case CUPPA_IMAGE:
                return new CuppaImageComparer(config, categoryFieldOverrides);

            case CHORD:
                return new ChordComparer(config, categoryFieldOverrides);

            case LILAC_QC:
                return new LilacQcComparer(config, categoryFieldOverrides);

            case LILAC_ALLELE:
                return new LilacAlleleComparer(config, categoryFieldOverrides);

            case GERMLINE_SV:
                return new GermlineSvComparer(config, categoryFieldOverrides);

            case PEACH:
                return new PeachComparer(config, categoryFieldOverrides);

            case VIRUS:
                return new VirusComparer(config, categoryFieldOverrides);

            case TUMOR_FLAGSTAT:
                return new FlagstatComparer(TUMOR_FLAGSTAT, config, categoryFieldOverrides);

            case GERMLINE_FLAGSTAT:
                return new FlagstatComparer(GERMLINE_FLAGSTAT, config, categoryFieldOverrides);

            case TUMOR_BAM_METRICS:
                return new BamMetricsComparer(TUMOR_BAM_METRICS, config, categoryFieldOverrides);

            case GERMLINE_BAM_METRICS:
                return new BamMetricsComparer(GERMLINE_BAM_METRICS, config, categoryFieldOverrides);

            case SNP_GENOTYPE:
                return new SnpGenotypeComparer(config, categoryFieldOverrides);

            case CDR3_SEQUENCE:
                return new CiderVdjComparer(config, categoryFieldOverrides);

            case CDR3_LOCUS_SUMMARY:
                return new Cdr3LocusSummaryComparer(config, categoryFieldOverrides);

            case TELOMERE_LENGTH:
                return new TealComparer(config, categoryFieldOverrides);

            case V_CHORD:
                return new VChordComparer(config, categoryFieldOverrides);

            case SIGS:
                return new SigsComparer(config, categoryFieldOverrides);

            case RNA_SUMMARY:
                return new RnaSummaryComparer(config, categoryFieldOverrides);

            case RNA_GENE_DATA:
                return new RnaGeneDataComparer(config, categoryFieldOverrides);

            case RNA_TRANSCRIPT_DATA:
                return new RnaTranscriptDataComparer(config, categoryFieldOverrides);

            case RNA_NOVEL_SPLICE_JUNCTION:
                return new RnaNovelSpliceJunctionComparer(config, categoryFieldOverrides);

            case RNA_FUSION:
                return new RnaFusionComparer(config, categoryFieldOverrides);

            default:
                return null;
        }
    }
}
