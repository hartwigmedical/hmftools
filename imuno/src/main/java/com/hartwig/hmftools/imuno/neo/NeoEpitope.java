package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.common.neo.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.adjustCodingBasesForStrand;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.calcNonMediatedDecayBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.calcStartCodonBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.checkTrimBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.getAminoAcids;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.TranscriptCodingType;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;

import htsjdk.samtools.Cigar;

public abstract class NeoEpitope
{
    // transcript context
    public final TranscriptData[] TransData; // only start is populated for same-gene NEs

    public final int[] Phases;
    public final int[] ExonRank;

    public TranscriptCodingType[] CodingType;
    public TranscriptRegionType[] RegionType;

    public final String[] RawCodingBases; // coding bases from before and after the mutation or fusion junction
    public final String[] CodingBases; // corrected for strand and stripped of any novel amino-acid bases
    public final int[] NovelBaseIndex; // index into the up and downstream coding bases
    public String NovelCodonBases;

    public String UpstreamAcids;
    public String DownstreamAcids;
    public String NovelAcid;
    public int NmdBasesMin;
    public int NmdBasesMax;
    public int StartCodonBasesMin;
    public int StartCodonBasesMax;
    public String WildtypeAcids;

    public boolean Valid;

    // data for RNA matching
    public String[] ExtCodingBases;
    public int[][] ExtPositions; // coding base up and down position
    public final Cigar[] ExtCigars;

    public NeoEpitope()
    {
        TransData = new TranscriptData[] {null, null};
        Phases = new int[] {-1, -1};
        ExonRank = new int[FS_PAIR];
        CodingType = new TranscriptCodingType[] {TranscriptCodingType.UNKNOWN, TranscriptCodingType.UNKNOWN};
        RegionType = new TranscriptRegionType[] {TranscriptRegionType.UNKNOWN, TranscriptRegionType.UNKNOWN};

        RawCodingBases = new String[] {"", ""};
        CodingBases = new String[] {"", ""};
        NovelBaseIndex = new int[] {0, 0};
        NovelCodonBases = "";
        UpstreamAcids = "";
        DownstreamAcids = "";
        NovelAcid = "";
        NmdBasesMin = 0;
        NmdBasesMax = 0;
        StartCodonBasesMin = 0;
        StartCodonBasesMax = 0;
        WildtypeAcids = "";
        Valid = true;

        ExtCodingBases = new String[] {"", ""};
        ExtCigars = new Cigar[FS_PAIR];
        ExtPositions = new int[FS_PAIR][SE_PAIR];
    }

    public byte orientation(int fs)
    {
        return (fs == FS_UP) == (TransData[fs].Strand == POS_STRAND) ? POS_ORIENT : NEG_ORIENT;
    }

    public byte strand(int stream) { return TransData[stream].Strand; }

    public static int getUpstreamOpenCodonBases(int phase)
    {
        // assumes a codon starts with phase 1
        return phase == PHASE_NONE ? 0 : phase;
    }

    public abstract int position(int stream);
    public abstract String chromosome(int stream);
    public abstract String geneName(int stream);
    public abstract void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData);
    public abstract void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids);
    public abstract NeoEpitopeType variantType();
    public abstract String variantInfo();
    public abstract double copyNumber();
    public abstract boolean phaseMatched();
    public abstract int unsplicedDistance();
    public abstract int skippedAcceptors();
    public abstract int skippedDonors();

    public void setCodingBases(final RefGenomeInterface refGenome, int reqAminoAcids)
    {
        extractCodingBases(refGenome, reqAminoAcids);

        adjustCodingBasesForStrand(this);

        boolean isPhased = phaseMatched();

        int novelUpstreamBases = NovelBaseIndex[FS_UP];
        int novelDownstreamBases = NovelBaseIndex[FS_DOWN];

        if(novelUpstreamBases < 0 || novelDownstreamBases < 0)
        {
            IM_LOGGER.error("ne({}) invalid upBases({} noveBases={}) or downBases({} noveBases={})",
                    this, CodingBases[FS_UP], novelUpstreamBases, CodingBases[FS_DOWN], novelDownstreamBases);
            return;
        }

        // if upstream ends on a phase other than 2, need to take the bases from the downstream gene to make a novel codon
        if(novelUpstreamBases > 0 || novelDownstreamBases > 0 || !isPhased)
        {
            String upstreamBases = CodingBases[FS_UP];
            String downstreamBases = CodingBases[FS_DOWN];

            if(novelUpstreamBases > upstreamBases.length() || novelDownstreamBases > downstreamBases.length())
            {
                IM_LOGGER.error("ne({}) invalid upBases({} noveBases={}) or downBases({} noveBases={})",
                        this, upstreamBases, novelUpstreamBases, downstreamBases, novelDownstreamBases);
                return;
            }

            // take the last 1 or 2 bases from the end of upstream gene's section
            NovelCodonBases = upstreamBases.substring(upstreamBases.length() - novelUpstreamBases);
            upstreamBases = upstreamBases.substring(0, upstreamBases.length() - novelUpstreamBases);

            if(isPhased)
            {
                NovelCodonBases += downstreamBases.substring(0, novelDownstreamBases);
                downstreamBases = downstreamBases.substring(novelDownstreamBases);
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

            CodingBases[FS_UP] = upstreamBases;
            CodingBases[FS_DOWN] = downstreamBases;
        }

        IM_LOGGER.trace("ne({}) upBases({}) novelCodon({}) downBases({}) downNmdBases({})",
                this, CodingBases[FS_UP], checkTrimBases(NovelCodonBases),
                checkTrimBases(CodingBases[FS_DOWN]), NmdBasesMin);
    }

    public void setNonsenseMediatedDecay()
    {
        NmdBasesMin = NmdBasesMax = calcNonMediatedDecayBases(this);

        // also record the bases to the start codon
        StartCodonBasesMin = StartCodonBasesMax = calcStartCodonBases(this);
    }

    public void setAminoAcids()
    {
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

    public abstract void setSkippedSpliceSites(final EnsemblDataCache geneTransCache);

    public boolean matchesAminoAcids(final NeoEpitope other)
    {
        return UpstreamAcids.equals(other.UpstreamAcids) && DownstreamAcids.equals(other.DownstreamAcids) && NovelAcid.equals(other.NovelAcid);
    }

    public String aminoAcidString() { return UpstreamAcids + NovelAcid + DownstreamAcids; }

    public NeoEpitopeFile toFile(
            final Set<String> upTransNames, final Set<String> downTransNames,
            final double[] tpmCancer, final double[] tpmCohort)
    {
        final StringJoiner upTransStr = new StringJoiner(";");
        final StringJoiner downTransStr = new StringJoiner(";");
        upTransNames.forEach(x -> upTransStr.add(x));
        downTransNames.forEach(x -> downTransStr.add(x));

        return new NeoEpitopeFile(
                variantType(), variantInfo(), copyNumber(),
                TransData[FS_UP].GeneId, TransData[FS_DOWN].GeneId, geneName(FS_UP), geneName(FS_DOWN),
                chromosome(FS_UP), chromosome(FS_DOWN), orientation(FS_UP), orientation(FS_DOWN),
                UpstreamAcids, DownstreamAcids, NovelAcid, NmdBasesMin, NmdBasesMax, StartCodonBasesMin, StartCodonBasesMax,
                unsplicedDistance(), skippedDonors(), skippedAcceptors(),
                upTransStr.toString(), downTransStr.toString(), WildtypeAcids,
                ExtPositions[FS_UP][SE_START], ExtPositions[FS_UP][SE_END], ExtCodingBases[FS_UP], ExtCigars[FS_UP].toString(),
                ExtPositions[FS_DOWN][SE_START], ExtPositions[FS_DOWN][SE_END], ExtCodingBases[FS_DOWN],
                ExtCigars[FS_DOWN] != null ? ExtCigars[FS_DOWN].toString() : "",
                tpmCancer[FS_UP], tpmCohort[FS_UP], tpmCancer[FS_DOWN], tpmCohort[FS_DOWN]);
    }
}
