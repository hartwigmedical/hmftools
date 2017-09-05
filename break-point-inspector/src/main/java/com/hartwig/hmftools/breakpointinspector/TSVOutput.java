package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

import org.apache.commons.lang3.ObjectUtils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

class TSVOutput {

    private static List<String> parseMantaPRSR(final Genotype genotype) {
        String pr = (String) genotype.getExtendedAttribute("PR", "0,0");
        String sr = (String) genotype.getExtendedAttribute("SR", "0,0");
        return Stream.concat(Arrays.stream(pr.split(",")), Arrays.stream(sr.split(","))).collect(Collectors.toList());
    }

    static void PrintHeaders() {
        final ArrayList<String> header =
                Lists.newArrayList("ID", "SVTYPE", "ORIENTATION", "MANTA_BP1", "MANTA_BP2", "MANTA_SVLEN", "MANTA_REF_PR_NORMAL",
                        "MANTA_REF_PR_SUPPORT", "MANTA_REF_SR_NORMAL", "MANTA_REF_SR_SUPPORT", "MANTA_TUMOR_PR_NORMAL",
                        "MANTA_TUMOR_PR_SUPPORT", "MANTA_TUMOR_SR_NORMAL", "MANTA_TUMOR_SR_SUPPORT", "MANTA_HOMSEQ", "MANTA_INSSEQ");
        header.addAll(prefixList(SampleStats.GetHeader(), "REF_"));
        header.addAll(prefixList(SampleStats.GetHeader(), "TUMOR_"));
        header.add("BPI_BP1");
        header.add("BPI_BP2");
        header.add("FILTER");
        header.add("AF_BP1");
        header.add("AF_BP2");
        header.add("EXTRA_UNCERTAINTY");
        System.out.println(String.join("\t", header));
    }

    static void PrintVariant(final VariantContext variant, final HMFVariantContext context, final StructuralVariantResult result) {
        final List<String> fields = Lists.newArrayList(variant.getID(), variant.getStructuralVariantType().toString(),
                HMFVariantType.getOrientation(context.Type), result.toString(), context.MantaBP1.toString(),
                variant.getAttributeAsString("SVLEN", ""));

        fields.addAll(parseMantaPRSR(variant.getGenotype(0)));
        fields.addAll(parseMantaPRSR(variant.getGenotype(1)));

        fields.add(variant.getAttributeAsString("HOMSEQ", ""));
        fields.add(variant.getAttributeAsString("SVINSSEQ", ""));

        fields.addAll(result.RefStats.GetData());
        fields.addAll(result.TumorStats.GetData());
        fields.add(ObjectUtils.firstNonNull(result.Breakpoints.getLeft(), "err").toString());
        fields.add(ObjectUtils.firstNonNull(result.Breakpoints.getRight(), "err").toString());
        fields.add(result.FilterString);

        fields.add(String.format("%.2f", result.AlleleFrequency.getLeft()));
        fields.add(String.format("%.2f", result.AlleleFrequency.getRight()));
        fields.add(Integer.toString(result.ExtraUncertainty));

        System.out.println(String.join("\t", fields));
    }
}
