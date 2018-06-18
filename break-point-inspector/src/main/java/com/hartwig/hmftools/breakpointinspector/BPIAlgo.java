package com.hartwig.hmftools.breakpointinspector;

import static java.util.Arrays.asList;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

final class BPIAlgo {

    private BPIAlgo() {
    }

    @NotNull
    static BPIAlgoOutput run(@NotNull final VCFFileReader vcfReader, @NotNull final SAMSequenceDictionary dictionary,
            @NotNull final Analysis analysis) throws IOException {
        final List<VariantContext> variants = Lists.newArrayList();
        final List<QueryInterval> combinedQueryIntervals = Lists.newArrayList();
        final List<String> tsv = Lists.newArrayList();

        final Map<String, VariantContext> variantPerIDMap = Maps.newHashMap();

        for (VariantContext variant : vcfReader) {
            variantPerIDMap.put(variant.getID(), variant);

            final VariantContext mateVariant = variant;
            if (variant.hasAttribute("MATEID")) {
                variant = variantPerIDMap.get(variant.getAttributeAsString("MATEID", ""));
                if (variant == null) {
                    continue;
                }
            }

            final EnrichedVariantContext enrichedVariant = VariantEnrichment.enrich(variant, mateVariant, dictionary);

            final StructuralVariantResult result = analysis.processStructuralVariant(enrichedVariant);
            combinedQueryIntervals.addAll(asList(result.QueryIntervals));

            tsv.add(TSVOutput.generateVariant(enrichedVariant, result));

            final BiConsumer<VariantContext, Boolean> vcfUpdater = (v, swap) -> {
                final Set<String> filters = v.getCommonInfo().getFiltersMaybeNull();
                if (filters != null) {
                    filters.clear();
                }
                // NERA: We will map BreakpointError to a flag
                if (result.Filters.contains(FilterType.BREAKPOINT_ERROR.toString())) {
                    v.getCommonInfo().putAttribute("BPI_AMBIGUOUS", true, true);
                } else {
                    v.getCommonInfo().addFilters(result.Filters);
                }
                if (result.Filters.isEmpty()) {
                    final List<Double> af = asList(result.AlleleFrequency.getLeft(), result.AlleleFrequency.getRight());
                    v.getCommonInfo().putAttribute(AlleleFrequency.VCF_AF_INFO_TAG, swap ? Lists.reverse(af) : af, true);
                }
                if (result.Breakpoints.getLeft() != null) {
                    v.getCommonInfo().putAttribute(swap ? "BPI_END" : "BPI_START", result.Breakpoints.getLeft().position(), true);
                }
                if (result.Breakpoints.getRight() != null) {
                    v.getCommonInfo().putAttribute(swap ? "BPI_START" : "BPI_END", result.Breakpoints.getRight().position(), true);
                }

                // NERA: Remove CIPOS / CIEND when we have an insert sequence
                if (!v.hasAttribute("IMPRECISE") && v.hasAttribute("SVINSSEQ")) {
                    v.getCommonInfo().removeAttribute("CIPOS");
                    v.getCommonInfo().removeAttribute("CIEND");
                }
                variants.add(v);
            };

            vcfUpdater.accept(variant, false);
            if (mateVariant != variant) {
                vcfUpdater.accept(mateVariant, true);
            }
        }

        final QueryInterval[] optimizedIntervals =
                QueryInterval.optimizeIntervals(combinedQueryIntervals.toArray(new QueryInterval[combinedQueryIntervals.size()]));

        return ImmutableBPIAlgoOutput.builder().optimizedIntervals(optimizedIntervals).variants(variants).tsvOutput(tsv).build();
    }
}
