package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.LINX_DATA;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class LinxSvData implements ComparableItem
{
    public final StructuralVariantData SvData;
    public final LinxSvAnnotation Annotation;
    public final LinxCluster Cluster;

    public LinxSvData(
            final StructuralVariantData svData, final LinxSvAnnotation annotation, final LinxCluster cluster)
    {
        SvData = svData;
        Cluster = cluster;
        Annotation = annotation;
    }

    @Override
    public Category category() { return LINX_DATA; }

    @Override
    public String key()
    {
        return String.format("%d_%s %s_%d - %s_%d",
                SvData.id(), SvData.type(), SvData.startChromosome(), SvData.startPosition(), SvData.endChromosome(), SvData.endPosition());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
//        values.add(String.format("Qual(%.0f)", Variant.qual()));
//        values.add(String.format("Tier(%s)", Variant.tier().toString()));
//        values.add(String.format("TotalReadCount(%d)", Variant.totalReadCount()));
//        values.add(String.format("AlleleReadCount(%d)", Variant.alleleReadCount()));
        return values;
    }

    @Override
    public boolean reportable() { return false; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final LinxSvData otherSv = (LinxSvData)other;

        if(otherSv.SvData.type() != SvData.type())
            return false;

        if(!otherSv.SvData.startChromosome().equals(SvData.startChromosome()) || !otherSv.SvData.endChromosome().equals(SvData.endChromosome()))
            return false;

        if(otherSv.SvData.startPosition() != SvData.startPosition() || otherSv.SvData.endPosition() != SvData.endPosition())
            return false;

        if(otherSv.SvData.startOrientation() != SvData.startOrientation() || otherSv.SvData.endOrientation() != SvData.endOrientation())
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel)
    {
        final LinxSvData otherSv = (LinxSvData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "resolvedType", Cluster.resolvedType(), otherSv.Cluster.resolvedType());

        if(matchLevel == REPORTABLE)
            return null;

        checkDiff(diffs, "clusterCount", Cluster.clusterCount(), otherSv.Cluster.clusterCount());

        double jcn = Annotation.junctionCopyNumberMin() + Annotation.junctionCopyNumberMax();
        double jcnOther = otherSv.Annotation.junctionCopyNumberMin() + otherSv.Annotation.junctionCopyNumberMax();

        checkDiff(diffs, "jcn", jcn, jcnOther);

        if(!genesEqual(Annotation.geneStart(), otherSv.Annotation.geneStart()) || !genesEqual(Annotation.geneEnd(), otherSv.Annotation.geneEnd()))
        {
            diffs.add(String.format("genes(%s-%s/%s-%s)",
                    Annotation.geneStart(), Annotation.geneEnd(), otherSv.Annotation.geneStart(), otherSv.Annotation.geneEnd()));
        }

        checkDiff(diffs, "foldback", Annotation.isFoldback(), otherSv.Annotation.isFoldback());

        return null;
    }

    private boolean genesEqual(final String genes1, final String genes2)
    {
        final List<String> geneNames1 = Arrays.stream(genes1.split(ITEM_DELIM, -1)).collect(Collectors.toList());
        final List<String> geneNames2 = Arrays.stream(genes2.split(ITEM_DELIM, -1)).collect(Collectors.toList());

        if(geneNames1.size() != geneNames2.size())
            return false;

        return !geneNames1.stream().noneMatch(x -> geneNames2.contains(x));
    }
}
