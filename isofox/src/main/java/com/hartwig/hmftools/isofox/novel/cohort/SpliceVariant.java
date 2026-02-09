package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.ACCEPTOR;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.OLD_FILE_DELIM;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
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
    public final String CodingImpact; // HGVS canonical coding impact
    public final String TriNucContext;
    public final int LocalPhaseSet;

    public List<String> SampleIds;

    public SpliceVariant(
            final String geneName, final String chromosome, int position, final VariantType type,
            final String ref, final String alt, final String codingEffect, final String codingImpact,
            final String triNucContext, int localPhaseSet)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        CodingImpact = codingImpact;
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
                items[fieldIndexMap.get("CodingImpact")],
                items[fieldIndexMap.get("TriNucContext")],
                fieldIndexMap.containsKey("LocalPhaseSet") ? Integer.parseInt(items[fieldIndexMap.get("LocalPhaseSet")]) : 0);

    }

    public static SpliceVariant fromCsv(final String[] items)
    {
        int index = 0;
        return new SpliceVariant(
                items[index++],
                items[index++],
                Integer.parseInt(items[index++]),
                VariantType.valueOf(items[index++]),
                items[index],
                items[index],
                items[index],
                items[index],
                items[index],
                Integer.parseInt(items[index++]));

    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(OLD_FILE_DELIM);
        sj.add(FLD_GENE_NAME);
        sj.add(FLD_CHROMOSOME);
        sj.add("Position");
        sj.add("Type");
        sj.add("Ref");
        sj.add("Alt");
        sj.add("CodingEffect");
        sj.add("CodingImpact");
        sj.add("TriNucContext");
        sj.add("LocalPhaseSet");
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(OLD_FILE_DELIM);
        sj.add(GeneName);
        sj.add(Chromosome);
        sj.add(String.valueOf(Position));
        sj.add(String.valueOf(Type));
        sj.add(Ref);
        sj.add(Alt);
        sj.add(CodingEffect);
        sj.add(CodingImpact);
        sj.add(TriNucContext);
        sj.add(String.valueOf(LocalPhaseSet));
        return sj.toString();
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

        int startBasesReq = acceptorDonorType == ACCEPTOR ? 10 : 1;
        int endBasesReq = acceptorDonorType == ACCEPTOR ? 1: 10;

        int varDiff = variantAlt.length() - variantRef.length();

        boolean varWithinContext = variantPos < splicePosition ?
                (splicePosition - variantPos <= startBasesReq) : ((variantPos - splicePosition <= endBasesReq));

        int preSpliceBaseDiff = 0;
        int postSpliceBasesRequired = endBasesReq;

        if(varWithinContext && varDiff != 0)
        {
            if(variantPos < splicePosition)
            {
                preSpliceBaseDiff = -varDiff;
            }
            else
            {
                if(varDiff < 0)
                {
                    postSpliceBasesRequired += abs(varDiff); // need extra bases to cover the DEL
                }
                else
                {
                    // for the INS case, depends on where the var is relative to the splice site
                    // endBasesReq = 1, means 2 bases would normally be retrieved starting at the splice site
                    // if var = splice site, and var is an insert such that alt length >= end offset then no bases are required
                    postSpliceBasesRequired = variantPos - splicePosition;

                    if(postSpliceBasesRequired < endBasesReq)
                        postSpliceBasesRequired += max(endBasesReq - postSpliceBasesRequired - varDiff, 0);
                }
            }
        }

        final String preSpliceBases = startBasesReq + preSpliceBaseDiff >= 1 ?
                refGenome.getBaseString(chromosome, splicePosition - startBasesReq - preSpliceBaseDiff, splicePosition - 1) : "";

        final String postSpliceBases = postSpliceBasesRequired >= 0 ?
                refGenome.getBaseString(chromosome, splicePosition, splicePosition + postSpliceBasesRequired) : "";

        if(!varWithinContext)
            return preSpliceBases + postSpliceBases;

        String baseContext = "";
        int varRelativePos = 0;

        try
        {
            if(variantPos < splicePosition)
            {
                varRelativePos = preSpliceBases.length() - (splicePosition - variantPos);
                baseContext = preSpliceBases.substring(0, varRelativePos);
                baseContext += variantAlt;
                baseContext += preSpliceBases.substring(varRelativePos + variantRef.length(), preSpliceBases.length());
                baseContext += postSpliceBases;
            }
            else
            {
                baseContext = preSpliceBases;
                varRelativePos = variantPos - splicePosition;
                ++endBasesReq; // to account for the splice base

                if(postSpliceBases.isEmpty())
                {
                    baseContext += variantAlt.length() > endBasesReq ? variantAlt.substring(0, endBasesReq) : variantAlt;
                }
                else
                {
                    baseContext += postSpliceBases.substring(0, varRelativePos);

                    if(varRelativePos + variantAlt.length() > endBasesReq)
                    {
                        baseContext += variantAlt.substring(0, endBasesReq - varRelativePos);
                    }
                    else
                    {
                        baseContext += variantAlt;
                        baseContext += postSpliceBases.substring(varRelativePos + variantRef.length(), postSpliceBases.length());
                    }
                }
            }
        }
        catch (Exception e)
        {
            ISF_LOGGER.error("base context error: {}", e.toString());

            ISF_LOGGER.error("var(chr={} pos={} ref={} alt={}) splice(pos={} ad={})",
                    chromosome, variantPos, variantRef, variantAlt, splicePosition, acceptorDonorType);

            ISF_LOGGER.error("baseContext({}) varRelPos({}) preSpliceBases({}) postSpliceBases({})",
                    baseContext, varRelativePos, preSpliceBases, postSpliceBases);
        }

        return baseContext;
    }

    public static String formKey(final String chromsome, final int position, final String ref, final String alt)
    {
        return String.format("%s-%d-%s-%s", chromsome, position, ref, alt);
    }

    public String key() { return formKey(Chromosome, Position, Ref, Alt); }

    public boolean matchesLocation(final String chromosome, final int position)
    {
        return Chromosome.equals(chromosome) && Position == position;
    }

    public boolean matches(final String chromosome, final int position, final String ref, final String alt)
    {
        return matchesLocation(chromosome, position) && Ref.equals(ref) && Alt.equals(alt);
    }

    public void addSampleId(final String sampleId)
    {
        if(SampleIds == null)
        {
            SampleIds = Lists.newArrayList();
        }

        if(!SampleIds.contains(sampleId))
            SampleIds.add(sampleId);
    }
}
