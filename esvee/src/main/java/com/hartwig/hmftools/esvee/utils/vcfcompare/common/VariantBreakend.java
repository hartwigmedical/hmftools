package com.hartwig.hmftools.esvee.utils.vcfcompare.common;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.EVENT_TYPE;
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

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.utils.vcfcompare.line.LineLinker;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantBreakend
{
    public final VariantContext Context;
    private final VariantAltInsertCoords AltCoords;

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
    public VariantBreakend LinkedLineBreakend;

    public VariantBreakend(final VariantContext context, SvCaller svCaller, VcfType sourceVcfType)
    {
        Context = context;

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

        SvType = svCaller == SvCaller.GRIDSS ?
                context.getAttributeAsString(EVENT_TYPE, "") :
                context.getAttributeAsString(SV_TYPE, "");

        Filters = Context.getFilters();

        SourceVcfType = sourceVcfType;

        MatchedBreakend = null;
        LinkedLineBreakend = null;
    }

    public boolean isSingle()
    {
        return OtherChromosome.isEmpty();
    }

    public boolean isTranslocation()
    {
        return !isSingle() && Chromosome.equals(OtherChromosome);
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

    public boolean isLineInsertionSite() { return LineLinker.isLineInsertionSite(this); }

    public boolean hasLineLink() { return LinkedLineBreakend != null; }

    @Deprecated
    public boolean hasSglLineLink() { return hasLineLink() && isSingle() && LinkedLineBreakend.isSingle(); }

    public String toString()
    {
        return String.format("coords(%s) cipos(%d,%d)", coordStr(), Cipos[0], Cipos[1]);
    }

    public String coordStr()
    {
        return String.format("%s:%d:%d", Context.getContig(), Position, Orientation);
    }

    public String otherCoordStr()
    {
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

    public String lineLinkCoordStr()
    {
        if(!hasLineLink())
            return coordStr();

        if(isLowerBreakend(LinkedLineBreakend))
            return coordStr() + "_" + LinkedLineBreakend.coordStr();
        else
            return LinkedLineBreakend.coordStr() + "_" + coordStr();
    }

    public String filtersStr()
    {
        return isPassVariant() ? PASS : String.join(",", Filters);
    }

    public String qualStr()
    {
        return String.format("%.0f", Context.getPhredScaledQual());
    }

    public String id() { return Context.getID(); }

    public String fragsStr(String sampleId){ return getExtendedAttributeAsString(sampleId, TOTAL_FRAGS); }

    public StructuralVariantType svType()
    {
        if(OtherChromosome.equals(""))
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

    public static Map<String,List<VariantBreakend>> loadVariants(final String vcfFile)
    {
        SV_LOGGER.info("Loading vcfFile({})", vcfFile);

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

        SV_LOGGER.debug("  Loaded {} SVs from {})",
                chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum(),
                vcfFile
        );

        return chrBreakendMap;
    }
}
