package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.compar.common.CategoryType.DISRUPTION;
import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.FULL_MATCH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_ERROR;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_OLD;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.io.File;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.chord.ChordComparer;
import com.hartwig.hmftools.compar.cider.Cdr3LocusSummaryComparer;
import com.hartwig.hmftools.compar.cider.CiderVdjComparer;
import com.hartwig.hmftools.compar.common.field.Field;
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

public class CommonUtils
{
    public static final String FLD_REPORTED = "Reported";
    public static final String FLD_QUAL = "Qual";
    public static final String FLD_CHROMOSOME_BAND = "ChromosomeBand";
    public static final String FLD_FILTER = "filter";

    public static List<ItemComparer> buildComparers(final ComparConfig config)
    {
        List<ItemComparer> comparers = Lists.newArrayList();

        // load in a predictable order irrespective of config
        for(CategoryType category : CategoryType.values())
        {
            if(!config.Categories.containsKey(category))
            {
                continue;
            }

            MatchLevel matchLevel = config.Categories.get(category);

            ItemComparer comparer = createComparer(category, config);

            if(matchLevel == REPORTABLE && !comparer.hasReportable())
            {
                continue;
            }

            comparers.add(comparer);
        }

        // link related or dependent comparers - could make this a virtual method too if becomes more common
        if(config.Categories.containsKey(FUSION))
        {
            FusionComparer fusionComparer = (FusionComparer)(comparers.stream()
                    .filter(x -> x.category() == FUSION).findFirst().orElse(null));

            DisruptionComparer disruptionComparer;
            if(config.Categories.containsKey(DISRUPTION))
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

    public static boolean processSample(final ItemComparer comparer, final ComparConfig config, final String sampleId,
            final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        final MatchLevel matchLevel = config.Categories.get(comparer.category());

        Map<SourceType,List<ComparableItem>> sourceItems = Maps.newHashMap();

        for(SourceData source : config.Sources)
        {
            String sourceSampleId = config.sourceSampleId(source.Type, sampleId);
            String sourceReferenceId = config.sourceReferenceId(source.Type, sampleId);
            List<ComparableItem> items = null;

            if(source.Database != null)
            {
                items = comparer.loadFromDb(sourceSampleId, source.Database, source.Type);
            }
            else
            {
                FileSources fileSources = FileSources.sampleInstance(source.Files, sourceSampleId, sourceReferenceId);

                items = comparer.loadFromFile(sourceSampleId, sourceReferenceId, fileSources);
            }

            if(items != null)
            {
                sourceItems.put(source.Type, items);
            }
        }

        if(sourceItems.containsKey(SourceType.OLD) && sourceItems.containsKey(SourceType.NEW))
        {
            // previously support comparisons for N sources but now can only be 2 as controlled by config
            CommonUtils.compareItems(
                    mismatches, matchLevel, fieldConfig, config.IncludeMatches,
                    sourceItems.get(SourceType.OLD), sourceItems.get(SourceType.NEW));

            return true;
        }

        InvalidDataItem invalidDataItem = new InvalidDataItem(comparer.category());

        if(!sourceItems.containsKey(SourceType.OLD) && !sourceItems.containsKey(SourceType.NEW))
        {
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_BOTH, Collections.EMPTY_LIST));
        }
        else if(!sourceItems.containsKey(SourceType.OLD))
        {
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_OLD, Collections.EMPTY_LIST));
        }
        else if(!sourceItems.containsKey(SourceType.NEW))
        {
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_NEW, Collections.EMPTY_LIST));
        }

        return false;
    }

    public static void compareItems(
            final List<Mismatch> mismatches, final MatchLevel matchLevel, final FieldConfig fieldConfig, final boolean includeMatches,
            final List<ComparableItem> items1, final List<ComparableItem> items2)
    {
        int index1 = 0;
        while(index1 < items1.size())
        {
            final ComparableItem item1 = items1.get(index1);

            boolean matched = false;

            int index2 = 0;
            while(index2 < items2.size())
            {
                final ComparableItem item2 = items2.get(index2);

                if(item1.matches(item2))
                {
                    items1.remove(index1);
                    items2.remove(index2);
                    matched = true;

                    // skip checking for diffs if the items are not reportable
                    boolean eitherReportable = item1.reportable() || item2.reportable();

                    if(matchLevel != REPORTABLE || eitherReportable)
                    {
                        Mismatch mismatch = item1.findMismatch(item2, matchLevel, fieldConfig, includeMatches);

                        if(mismatch != null)
                        {
                            mismatches.add(mismatch);
                        }
                    }

                    break;
                }
                else
                {
                    ++index2;
                }
            }

            if(!matched)
            {
                ++index1;
            }
        }

        if(items1.isEmpty() && items2.isEmpty())
        {
            return;
        }

        List<String> emptyDiffs = Lists.newArrayList();

        items1.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x, null, OLD_ONLY, emptyDiffs)));

        items2.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(null, x, NEW_ONLY, emptyDiffs)));
    }

    public static List<String> findDiffs(final ComparableItem oldItem, final ComparableItem newItem, List<Field> fields)
    {
        final List<String> diffs = Lists.newArrayList();
        for(Field field : fields)
        {
            if(field.isCompared())
            {
                diffs.addAll(field.determineDiffs(oldItem, newItem));
            }
        }
        return diffs;
    }

    public static BasePosition determineComparisonGenomePosition(
            final String chromosome, final int position, final SourceType sourceType,
            final boolean requiresLiftover, final GenomeLiftoverCache liftoverCache)
    {
        if(requiresLiftover && liftoverCache.hasMappings())
        {
            RefGenomeVersion destVersion = sourceType == SourceType.OLD ? V38 : V37;
            int newPosition = liftoverCache.convertPosition(chromosome, position, destVersion);

            if(newPosition != UNMAPPED_POSITION)
            {
                return new BasePosition(destVersion.versionedChromosome(chromosome), newPosition);
            }
        }

        return new BasePosition(chromosome, position);
    }

    public static String determineComparisonChromosome(final String chromosome, final boolean requiresLiftover)
    {
        if(requiresLiftover)
        {
            return HumanChromosome.fromString(chromosome).name().substring(1);
        }
        else
        {
            return chromosome;
        }
    }

    public static boolean fileExists(final String filename)
    {
        return Files.exists(new File(filename).toPath());
    }

    public static Mismatch createMismatchFromDiffs(
            final ComparableItem refItem, final ComparableItem newItem, final List<String> diffs,
            final MatchLevel matchLevel, final boolean includeMatches)
    {
        if(diffs.isEmpty() && !includeMatches)
        {
            return null;
        }

        boolean refCountsAsCalled = countsAsCalled(refItem, matchLevel);
        boolean newCountsAsCalled = countsAsCalled(newItem, matchLevel);
        if(!refCountsAsCalled && !newCountsAsCalled)
        {
            // ignore unimportant differences
            return null;
        }
        else if(refCountsAsCalled && newCountsAsCalled && diffs.isEmpty() && !includeMatches)
        {
            // ignore perfect matches when not including matches
            return null;
        }
        else
        {
            MismatchType mismatchType;
            if(refCountsAsCalled && !newCountsAsCalled)
            {
                mismatchType = OLD_ONLY;
            }
            else if(!refCountsAsCalled && newCountsAsCalled)
            {
                mismatchType = NEW_ONLY;
            }
            else if(refCountsAsCalled && newCountsAsCalled && !diffs.isEmpty())
            {
                mismatchType = VALUE;
            }
            else if(refCountsAsCalled && newCountsAsCalled && includeMatches)
            {
                mismatchType = FULL_MATCH;
            }
            else
            {
                // should be impossible due to earlier filters
                mismatchType = INVALID_ERROR;
            }
            return new Mismatch(refItem, newItem, mismatchType, diffs);
        }
    }

    public static boolean countsAsCalled(final ComparableItem item, final MatchLevel matchLevel)
    {
        if(matchLevel == REPORTABLE)
        {
            return item.reportable();
        }
        else
        {
            return item.isPass();
        }
    }

    public static FieldConfig initialiseFieldConfig(ComparConfig config)
    {
        FieldConfig fieldConfig = new FieldConfig();
        for(CategoryType category : CategoryType.values())
        {
            MatchLevel matchLevel = config.Categories.get(category);

            ItemComparer comparer = createComparer(category, config);

            fieldConfig.registerFields(comparer, matchLevel);
        }

        if(config.FieldOverrides != null)
        {
            fieldConfig.applyOverrides(config.FieldOverrides, config.StrictFieldConfig);
        }

        return fieldConfig;
    }
}
