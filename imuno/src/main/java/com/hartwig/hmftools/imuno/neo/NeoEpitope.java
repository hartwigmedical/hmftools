package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.imuno.neo.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.adjustCodingBasesForStrand;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.calcNonMediatedDecayBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.checkTrimBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.getAminoAcids;

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
        return (fs == FS_UP) == (TransData[fs].Strand == POS_STRAND) ? POS_ORIENT : NEG_ORIENT;
    }

    public byte strand(int stream) { return TransData[stream].Strand; }

    public boolean phaseMatched()
    {
        if(RegionType[FS_UP] == EXONIC && RegionType[FS_DOWN] == EXONIC)
            return tickPhaseForward(Phases[FS_UP]) == Phases[FS_DOWN];
        else
            return Phases[FS_UP] == Phases[FS_DOWN];
    }

    public static int getUpstreamOpenCodonBases(int phase)
    {
        // assumes a codon starts with phase 1
        return phase == PHASE_NONE ? 0 : phase;

        /*
        if(phase == 0)
            return 1;
        else if(phase == 1)
            return 2;
        else
            return 0;
        */
    }

    public int getUpstreamOpenCodonBases() { return getUpstreamOpenCodonBases(Phases[FS_UP]); }

    public int getDownstreamPhaseOffset()
    {
        int upOpenBases = getUpstreamOpenCodonBases();
        return upOpenBases == 0 ? 0 : 3 - upOpenBases;
    }

    public abstract int position(int stream);
    public abstract String chromosome(int stream);
    public abstract String geneName(int stream);
    public abstract void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData);
    public abstract void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids);
    public abstract String variantType();
    public abstract String variantInfo();
    public abstract double copyNumber();

    public void setCodingBases(final RefGenomeInterface refGenome, int reqAminoAcids)
    {
        extractCodingBases(refGenome, reqAminoAcids);

        adjustCodingBasesForStrand(this);

        DownstreamNmdBases = calcNonMediatedDecayBases(this, FS_DOWN);

        boolean isPhased = phaseMatched();
        int upPhaseOffset = getUpstreamOpenCodonBases();
        int downPhaseOffset = getDownstreamPhaseOffset();

        // if upstream ends on a phase other than 2, need to take the bases from the downstream gene to make a novel codon
        if(upPhaseOffset > 0 || downPhaseOffset > 0 || !isPhased)
        {
            String upstreamBases = CodingBases[FS_UP];
            String downstreamBases = CodingBases[FS_DOWN];

            if(upPhaseOffset > upstreamBases.length() || downPhaseOffset > downstreamBases.length())
            {
                IM_LOGGER.error("ne({}) invalid upBases({} phaseOffset={}) or downBases({} phaseOffset={})",
                        this, upstreamBases, upPhaseOffset, downstreamBases, downPhaseOffset);
                return;
            }

            // take the last 1 or 2 bases from the end of upstream gene's section
            NovelCodonBases = upstreamBases.substring(upstreamBases.length() - upPhaseOffset);
            upstreamBases = upstreamBases.substring(0, upstreamBases.length() - upPhaseOffset);

            if(isPhased)
            {
                NovelCodonBases += downstreamBases.substring(0, downPhaseOffset);
                downstreamBases = downstreamBases.substring(downPhaseOffset);
            }
            else
            {
                NovelCodonBases += downstreamBases;
                downstreamBases = "";
            }

            // check superflous bases from longer inserts
            int requiredLength = reqAminoAcids * 3;

            if(upstreamBases.length() > requiredLength)
                upstreamBases = upstreamBases.substring(upstreamBases.length() - requiredLength);

            /*
            if(downstreamBases.length() > requiredLength)
                downstreamBases = downstreamBases.substring(0, requiredLength);

            if(NovelCodonBases.length() > requiredLength + 3)
                NovelCodonBases = NovelCodonBases.substring(0, requiredLength + 3);
            */

            CodingBases[FS_UP] = upstreamBases;
            CodingBases[FS_DOWN] = downstreamBases;
        }

        IM_LOGGER.trace("ne({}) upBases({}) novelCodon({}) downBases({}) downNmdBases({})",
                this, CodingBases[FS_UP], checkTrimBases(NovelCodonBases),
                checkTrimBases(CodingBases[FS_DOWN]), DownstreamNmdBases);
    }

    public void setAminoAcids()
    {
        boolean isPhased = phaseMatched();
        UpstreamAcids = getAminoAcids(CodingBases[FS_UP], false);
        NovelAcid = getAminoAcids(NovelCodonBases, true);
        DownstreamAcids = getAminoAcids(CodingBases[FS_DOWN], true);

        // truncate the novel bases and AAs if a stop codon is encountered
        if(NovelAcid.endsWith(STOP_SYMBOL))
        {
            int novelBases = NovelAcid.length() * 3;
            NovelCodonBases = NovelCodonBases.substring(0, novelBases);
        }

        IM_LOGGER.trace("ne({}) upAA({}) novel({}) downAA({})",
                this, UpstreamAcids, checkTrimBases(NovelAcid), checkTrimBases(DownstreamAcids));

        if(NovelAcid.equals(STOP_SYMBOL))
            DownstreamAcids = "";
    }

    public boolean matchesAminoAcids(final NeoEpitope other)
    {
        return UpstreamAcids.equals(other.UpstreamAcids) && DownstreamAcids.equals(other.DownstreamAcids) && NovelAcid.equals(other.NovelAcid);
    }
}
