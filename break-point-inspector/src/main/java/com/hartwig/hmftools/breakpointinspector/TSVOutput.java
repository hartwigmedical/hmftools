package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;
import com.hartwig.hmftools.breakpointinspector.datamodel.SampleStats;
import com.hartwig.hmftools.breakpointinspector.datamodel.StructuralVariantResult;

import org.apache.commons.lang3.ObjectUtils;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;

final class TSVOutput {

    private TSVOutput() {
    }

    @NotNull
    static String generateHeaders() {
        final ArrayList<String> header = Lists.newArrayList("ID",
                "SVTYPE",
                "ORIENTATION",
                "MANTA_BP1",
                "MANTA_BP2",
                "MANTA_SVLEN",
                "MANTA_REF_PR_NORMAL",
                "MANTA_REF_PR_SUPPORT",
                "MANTA_REF_SR_NORMAL",
                "MANTA_REF_SR_SUPPORT",
                "MANTA_TUMOR_PR_NORMAL",
                "MANTA_TUMOR_PR_SUPPORT",
                "MANTA_TUMOR_SR_NORMAL",
                "MANTA_TUMOR_SR_SUPPORT",
                "MANTA_HOMSEQ",
                "MANTA_INSSEQ");
        header.addAll(prefixList(SampleStats.header(), "REF_"));
        header.addAll(prefixList(SampleStats.header(), "TUMOR_"));
        header.add("BPI_BP1");
        header.add("BPI_BP2");
        header.add("FILTER");
        header.add("AF_BP1");
        header.add("AF_BP2");
        return String.join("\t", header);
    }

    @NotNull
    static String generateVariant(@NotNull EnrichedVariantContext context, @NotNull StructuralVariantResult result) {
        final List<String> fields = Lists.newArrayList(context.variant().getID(),
                context.variant().getStructuralVariantType().toString(),
                context.type().orientationString(),
                context.locationBP1().toString(),
                context.locationBP2().toString(),
                context.variant().getAttributeAsString("SVLEN", ""));

        fields.addAll(parseMantaPRSR(context.variant().getGenotype(0)));
        fields.addAll(parseMantaPRSR(context.variant().getGenotype(1)));

        fields.add(context.variant().getAttributeAsString("HOMSEQ", ""));
        fields.add(context.variant().getAttributeAsString("SVINSSEQ", ""));

        fields.addAll(result.refStats.data());
        fields.addAll(result.tumorStats.data());
        fields.add(ObjectUtils.firstNonNull(result.breakpoints.getLeft(), "err").toString());
        fields.add(ObjectUtils.firstNonNull(result.breakpoints.getRight(), "err").toString());
        fields.add(result.filterString);

        fields.add(String.format("%.2f", result.alleleFrequency.getLeft()));
        fields.add(String.format("%.2f", result.alleleFrequency.getRight()));

        return String.join("\t", fields);
    }

    @NotNull
    private static List<String> parseMantaPRSR(@NotNull Genotype genotype) {
        String pr = (String) genotype.getExtendedAttribute("PR", "0,0");
        String sr = (String) genotype.getExtendedAttribute("SR", "0,0");
        return Stream.concat(Arrays.stream(pr.split(",")), Arrays.stream(sr.split(","))).collect(Collectors.toList());
    }
}
