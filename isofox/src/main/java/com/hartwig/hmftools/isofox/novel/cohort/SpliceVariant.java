package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.ACCEPTOR;

import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.VariantType;

public class SpliceVariant
{
    public final String GeneName;
    public final String Chromosome;

    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;
    public final String CodingEffect;
    public final String HgvsCodingImpact;
    public final String TriNucContext;
    public final int LocalPhaseSet;

    public SpliceVariant(
            final String geneName, final String chromosome, int position, final VariantType type,
            final String ref, final String alt, final String codingEffect, final String hgvsCodingImpact,
            final String triNucContext, int localPhaseSet)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        HgvsCodingImpact = hgvsCodingImpact;
        TriNucContext = triNucContext;
        CodingEffect = codingEffect;
        LocalPhaseSet = localPhaseSet;
    }

    public static SpliceVariant fromCsv(final String[] items, final Map<String,Integer> fieldIndexMap)
    {
        return new SpliceVariant(
                items[fieldIndexMap.get("GeneName")],
                items[fieldIndexMap.get("Chromosome")],
                Integer.parseInt(items[fieldIndexMap.get("Position")]),
                VariantType.valueOf(items[fieldIndexMap.get("Type")]),
                items[fieldIndexMap.get("Ref")],
                items[fieldIndexMap.get("Alt")],
                items[fieldIndexMap.get("CodingEffect")],
                items[fieldIndexMap.get("HgvsImpact")],
                items[fieldIndexMap.get("TriNucContext")],
                fieldIndexMap.containsKey("LocalPhaseSet") ? Integer.parseInt(items[fieldIndexMap.get("LocalPhaseSet")]) : 0);

    }

    public String toString()
    {
        return String.format("gene(%s) loc(%s:%d) type(%s)", GeneName, Chromosome, Position, Type);
    }

    public static String getBaseContext(
            final String chromosome, int variantPos, final String variantRef, final String variantAlt,
            int splicePosition, final AcceptorDonorType acceptorDonorType, final RefGenomeInterface refGenome)
    {
        // first establish the bases around the splice junction, taking into account the ref to alt change

        int startOffset = acceptorDonorType == ACCEPTOR ? 10 : 1;
        int endOffset = acceptorDonorType == ACCEPTOR ? 1: 10;

        int varDiff = variantAlt.length() - variantRef.length();

        int preSpliceBaseDiff = 0;
        int postSpliceBaseDiff = 0;

        boolean varWithinContext = variantPos < splicePosition ?
                (splicePosition - variantPos <= startOffset) : ((variantPos - splicePosition <= endOffset));

        if(varWithinContext && varDiff != 0)
        {
            // bases added from an insert
            if(variantPos < splicePosition)
                preSpliceBaseDiff = -varDiff;
            else
                postSpliceBaseDiff = -varDiff;
        }

        final String preSpliceBases = startOffset + preSpliceBaseDiff >= 1 ?
                refGenome.getBaseString(chromosome, splicePosition - startOffset - preSpliceBaseDiff, splicePosition - 1) : "";

        final String postSpliceBases = endOffset + postSpliceBaseDiff > 0 ?
                refGenome.getBaseString(chromosome, splicePosition, splicePosition + endOffset + postSpliceBaseDiff) : "";

        String baseContext = "";

        if(!varWithinContext)
            return preSpliceBases + postSpliceBases;

        if(variantPos < splicePosition)
        {
            int varRelativePos = preSpliceBases.length() - (splicePosition - variantPos);
            baseContext = preSpliceBases.substring(0, varRelativePos);
            baseContext += variantAlt;
            baseContext += preSpliceBases.substring(varRelativePos + variantRef.length(), preSpliceBases.length());
            baseContext += postSpliceBases;
        }
        else
        {
            baseContext = preSpliceBases;
            int varRelativePos = variantPos - splicePosition;
            baseContext += postSpliceBases.substring(0, varRelativePos);
            baseContext += variantAlt;
            baseContext += postSpliceBases.substring(varRelativePos + variantRef.length(), postSpliceBases.length());
        }

        return baseContext;
    }
}
