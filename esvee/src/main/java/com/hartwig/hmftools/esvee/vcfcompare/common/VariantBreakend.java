package com.hartwig.hmftools.esvee.vcfcompare.common;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.EVENT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.EVENT_TYPE;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.PAR_ID;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SGL_FRAG_COUNT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.formSvType;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.vcfcompare.line.LineLink;
import com.hartwig.hmftools.esvee.vcfcompare.line.LineLinker;

import org.apache.commons.lang3.NotImplementedException;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantBreakend
{
    public final SvCaller mSvCaller;

    public final VariantContext Context;
    private final VariantAltInsertCoords AltCoords;

    public final String Id;

    public final String Chromosome;
    public final int Position;
    public final byte Orientation;

    public String OtherChromosome;
    public int OtherPosition;
    public byte OtherOrientation;

    public final int[] Cipos;
    public final int[] Ihompos;
    public final String Homseq;
    public final String InsertSequence;

    public final String SvType;
    public final Set<String> Filters;
    public final VcfType SourceVcfType;

    public VariantBreakend MatchedBreakend;
    public LineLink LinkedLineBreakends;
    public LineLink InferredLinkedLineBreakends;

    public VariantBreakend(final VariantContext context, SvCaller svCaller, VcfType sourceVcfType)
    {
        mSvCaller = svCaller;

        Context = context;

        Id = Context.getID();

        String alt = context.getAlternateAllele(0).getDisplayString();
        AltCoords = VariantAltInsertCoords.fromRefAlt(alt, alt.substring(0, 1));

        Chromosome = Context.getContig();
        Position = Context.getStart();
        Orientation = AltCoords.Orient.asByte();

        OtherChromosome = AltCoords.OtherChromsome;
        OtherPosition = AltCoords.OtherPosition;
        OtherOrientation = AltCoords.OtherOrient == null ? 0 : AltCoords.OtherOrient.asByte();

        Cipos = getPositionOffsets(CIPOS);
        Ihompos = getPositionOffsets(IHOMPOS);
        Homseq = context.getAttributeAsString(HOMSEQ, "");
        InsertSequence = AltCoords.InsertSequence;

        SvType = mSvCaller == SvCaller.GRIDSS ?
                context.getAttributeAsString(EVENT_TYPE, "") :
                context.getAttributeAsString(SV_TYPE, "");

        Filters = Context.getFilters();

        SourceVcfType = sourceVcfType;

        MatchedBreakend = null;
        LinkedLineBreakends = null;
        InferredLinkedLineBreakends = null;
    }

    public boolean isSingle()
    {
        return OtherChromosome.isEmpty();
    }

    public boolean isTranslocation()
    {
        return !isSingle() && !Chromosome.equals(OtherChromosome);
    }

    private boolean isEnd()
    {
        if(OtherChromosome.isEmpty())
        {
            return false;
        }

        if(Context.getContig().equals(OtherChromosome))
        {
            return Position > OtherPosition;
        }

        return HumanChromosome.lowerChromosome(OtherChromosome, Chromosome);
    }

    public boolean isInverted()
    {
        return Orientation == OtherOrientation;
    }

    private int[] getPositionOffsets(String vcfTag)
    {
        List<Integer> offsetsList = Context.getAttributeAsIntList(vcfTag, 0);
        return offsetsList.size() == 2 ?
                new int[] { offsetsList.get(0), offsetsList.get(1) } :
                new int[] { 0, 0 };
    }

    public int minPosition()
    {
        return Position + Cipos[0];
    }

    public int maxPosition()
    {
        return Position + Cipos[1];
    }

    public int otherMinPosition()
    {
        return isInverted() ?
                OtherPosition - Cipos[0] :
                OtherPosition + Cipos[0];
    }

    public int otherMaxPosition()
    {
        return isInverted() ?
                OtherPosition - Cipos[1] :
                OtherPosition + Cipos[1];
    }

    public String getExtendedAttributeAsString(String id, String key)
    {
        Object value = Context.getGenotype(id).getExtendedAttribute(key);
        return value != null ? value.toString() : "";
    }

    public int getExtendedAttributeAsInt(String id, String key)
    {
        return getGenotypeAttributeAsInt(Context.getGenotype(id), key, 0);
    }

    public double getExtendedAttributeAsDouble(String id, String key)
    {
        return getGenotypeAttributeAsDouble(Context.getGenotype(id), key, 0);
    }

    public boolean hasMatchedBreakend() { return MatchedBreakend != null; }

    public boolean hasPolyATail() { return LineLinker.hasPolyATail(this); }

    public boolean hasLineLink() { return LinkedLineBreakends != null; }

    public boolean hasInferredLineLink() { return InferredLinkedLineBreakends != null; }

    public boolean hasLineInfoFlag() { return Context.getAttributeAsBoolean(LINE_SITE, false); }

    public String toString()
    {
        return String.format("id(%s) coords(%s) cipos(%d,%d)", Id, coordStr(), Cipos[0], Cipos[1]);
    }

    public String coordStr()
    {
        return String.format("%s:%d:%d", Context.getContig(), Position, Orientation);
    }

    public String otherCoordStr()
    {
        if(isSingle())
            return "";

        return String.format("%s:%d:%d", OtherChromosome, OtherPosition, OtherOrientation);
    }

    public String svCoordStr()
    {
        if(isSingle())
        {
            return coordStr();
        }

        if(isEnd())
        {
            return otherCoordStr() + "_" + coordStr();
        }

        return coordStr() + "_" + otherCoordStr();
    }

    private boolean isLowerBreakend(VariantBreakend otherBreakend)
    {
        if(otherBreakend.Chromosome.isEmpty())
        {
            return true;
        }

        if(Context.getContig().equals(otherBreakend.Chromosome))
        {
            return Position < otherBreakend.Position;
        }

        return HumanChromosome.lowerChromosome(Chromosome, otherBreakend.Chromosome);
    }

    public String filtersStr()
    {
        return isPassVariant() ? PASS : String.join(",", Filters);
    }

    public double qual()
    {
        return Context.getPhredScaledQual();
    }

    public String qualStr()
    {
        return String.format("%.0f", qual());
    }

    public String fragsStr(String sampleId)
    {
        String fragsInfoTag;

        if(mSvCaller == SvCaller.GRIDSS)
        {
            fragsInfoTag = isSingle() ? SGL_FRAG_COUNT : TOTAL_FRAGS;
        }
        else
        {
            fragsInfoTag = TOTAL_FRAGS;
        }

        return getExtendedAttributeAsString(sampleId, fragsInfoTag);
    }

    public StructuralVariantType svType()
    {
        if(isSingle())
        {
            return SGL;
        }

        return formSvType(
                Chromosome, OtherChromosome,
                Position, OtherPosition,
                AltCoords.Orient, AltCoords.OtherOrient,
                InsertSequence.isEmpty()
        );
    }

    public boolean isPassVariant()
    {
        return Filters.isEmpty();
    }

    public String mateId()
    {
        return mSvCaller == SvCaller.GRIDSS ?
                Context.getAttributeAsString(PAR_ID, "") :
                Context.getAttributeAsString(MATE_ID, "");
    }

    public String eventId()
    {
        // This id links 2 breakends together into 1 structural event

        if(mSvCaller == SvCaller.GRIDSS)
        {
            return Context.getAttributeAsString(EVENT, "");
        }

        if(mSvCaller == SvCaller.TRUTH)
            return Context.getAttributeAsString("SVID", "");

        if(mSvCaller == SvCaller.ESVEE)
        {
            if(isSingle())
                return Id;

            Set<String> mateIds = Sets.newHashSet(Id, mateId());
            return String.join(",",mateIds);
        }

        throw new NotImplementedException("eventId() not implemented for sv caller: " + mSvCaller);
    }

    public static Map<String,List<VariantBreakend>> loadVariants(final String vcfFile)
    {
        SV_LOGGER.info("Loading VCF file: {}", vcfFile);

        Map<String,List<VariantBreakend>> chrBreakendMap = new HashMap<>();

        VcfFileReader reader = new VcfFileReader(vcfFile);

        String currentChr = "";
        List<VariantBreakend> breakends = null;

        SvCaller svCaller = SvCaller.fromVcfPath(vcfFile);
        VcfType sourceVcfType = VcfType.fromVcfPath(vcfFile);

        for(VariantContext variantContext : reader.iterator())
        {
            String chromosome = variantContext.getContig();

            if(!currentChr.equals(chromosome))
            {
                currentChr = chromosome;
                breakends = new ArrayList<>();
                chrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(new VariantBreakend(variantContext, svCaller, sourceVcfType));
        }

        SV_LOGGER.debug("Loaded {} structural variants", chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum());

        return chrBreakendMap;
    }
}
