package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.AMINO_ACID_REF_COUNT;
import static com.hartwig.hmftools.svtools.neo.NeoEpitopeAnnotator.IM_LOGGER;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.adjustCodingBasesForStrand;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.calcNonMediatedDecayBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.checkTrimBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getAminoAcids;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.TranscriptCodingType;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public abstract class NeoEpitope
{
    // transcript context
    public final TranscriptData[] TransData; // only start is populated for same-gene NEs

    public final int[] Phases;
    public final int[] ExonRank;

    public TranscriptCodingType[] CodingType;
    public TranscriptRegionType[] RegionType;

    public final String[] CodingBases; // coding bases from before and after the mutation or fusion junction
    public String NovelCodonBases;

    public String UpstreamAcids;
    public String DownstreamAcids;
    public String NovelAcid;
    public int DownstreamNmdBases;

    public NeoEpitope()
    {
        TransData = new TranscriptData[] {null, null};
        Phases = new int[] {-1, -1};
        ExonRank = new int[FS_PAIR];
        CodingType = new TranscriptCodingType[] {TranscriptCodingType.UNKNOWN, TranscriptCodingType.UNKNOWN};
        RegionType = new TranscriptRegionType[] {TranscriptRegionType.UNKNOWN, TranscriptRegionType.UNKNOWN};

        CodingBases = new String[] {"", ""};
        NovelCodonBases = "";
        UpstreamAcids = "";
        DownstreamAcids = "";
        NovelAcid = "";
        DownstreamNmdBases = 0;
    }

    public byte orientation(int fs)
    {
        return (fs == FS_UPSTREAM) == (TransData[fs].Strand == POS_STRAND) ? POS_ORIENT : NEG_ORIENT;
    }

    public byte strand(int stream) { return TransData[stream].Strand; }

    public boolean phaseMatched()
    {
        if(RegionType[FS_UPSTREAM] == EXONIC && RegionType[FS_DOWNSTREAM] == EXONIC)
            return ((Phases[FS_UPSTREAM] + 1) % 3) == (Phases[FS_DOWNSTREAM] % 3);
        else
            return Phases[FS_UPSTREAM] == Phases[FS_DOWNSTREAM];
    }

    public int getDownstreamPhaseOffset()
    {
        return Phases[FS_UPSTREAM] == 0 || !phaseMatched() ? 0 : 3 - Phases[FS_UPSTREAM];
    }

    public abstract int position(int stream);
    public abstract String chromosome(int stream);
    public abstract void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData);
    public abstract void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids);

    public void setCodingBases(final RefGenomeInterface refGenome)
    {
        boolean isPhased = phaseMatched();

        int upstreamPhaseOffset = Phases[FS_UPSTREAM];
        int downstreamPhaseOffset = getDownstreamPhaseOffset();

        IM_LOGGER.debug("ne({}) phased({}) upPhaseOffset({}) downPhaseOffset({})",
                toString(), isPhased, upstreamPhaseOffset, downstreamPhaseOffset);

        extractCodingBases(refGenome, AMINO_ACID_REF_COUNT);

        adjustCodingBasesForStrand(this);

        DownstreamNmdBases = calcNonMediatedDecayBases(this, FS_DOWNSTREAM);

        String upstreamBases = CodingBases[FS_UPSTREAM];
        String downstreamBases = CodingBases[FS_DOWNSTREAM];

        if(upstreamPhaseOffset > upstreamBases.length() || downstreamPhaseOffset > downstreamBases.length())
        {
            IM_LOGGER.error("ne({}) invalid upBases({} phaseOffset={}) or downBases({} phaseOffset={})",
                    this, upstreamBases, upstreamPhaseOffset, downstreamBases, downstreamPhaseOffset);
            return;
        }

        // if upstream ends on a phase other than 0, need to take the bases from the downstream gene to make a novel codon
        if(upstreamPhaseOffset > 0)
        {
            // take the last 1 or 2 bases from the end of upstream gene's section
            NovelCodonBases = upstreamBases.substring(upstreamBases.length() - upstreamPhaseOffset);
            upstreamBases = upstreamBases.substring(0, upstreamBases.length() - upstreamPhaseOffset);
        }

        if(isPhased)
        {
            NovelCodonBases += downstreamBases.substring(0, downstreamPhaseOffset);
            downstreamBases = downstreamBases.substring(downstreamPhaseOffset);
        }
        else
        {
            NovelCodonBases += downstreamBases;
            downstreamBases = "";
        }

        CodingBases[FS_UPSTREAM] = upstreamBases;
        CodingBases[FS_DOWNSTREAM] = downstreamBases;

        IM_LOGGER.debug("ne({}) upBases({}) novelCodon({}) downBases({}) downNmdBases({})",
                this, upstreamBases, checkTrimBases(NovelCodonBases), checkTrimBases(downstreamBases), DownstreamNmdBases);
    }

    public void setAminoAcids()
    {
        boolean isPhased = phaseMatched();
        UpstreamAcids = getAminoAcids(CodingBases[FS_UPSTREAM], false);
        NovelAcid = getAminoAcids(NovelCodonBases, !isPhased);
        DownstreamAcids = getAminoAcids(CodingBases[FS_DOWNSTREAM], !isPhased);

        IM_LOGGER.debug("ne({}) upAA({}) novel({}) downAA({})",
                this, UpstreamAcids, checkTrimBases(NovelAcid), checkTrimBases(DownstreamAcids));

        if(NovelAcid.equals(STOP_SYMBOL))
            DownstreamAcids = "";
    }
}
