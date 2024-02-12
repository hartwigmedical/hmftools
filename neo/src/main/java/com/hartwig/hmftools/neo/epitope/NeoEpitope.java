package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.common.codon.AminoAcidRna.STOP_SYMBOL;
import static com.hartwig.hmftools.neo.NeoCommon.transcriptsToStr;

import java.util.Set;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
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
    public int CodingBasesLengthMin;
    public int CodingBasesLengthMax;
    public String UpstreamWildTypeAcids; // wildtype AAs starting from the upstream AAs

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
        CodingBasesLengthMin = 0;
        CodingBasesLengthMax = 0;
        UpstreamWildTypeAcids = "";
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
    public abstract double variantCopyNumber();
    public abstract double copyNumber();
    public abstract double subclonalLikelihood();
    public abstract boolean phaseMatched();
    public abstract int unsplicedDistance();
    public abstract int skippedAcceptors();
    public abstract int skippedDonors();
    public abstract String wildtypeAcids();

    public void setCodingBases(final RefGenomeInterface refGenome, int reqAminoAcids)
    {
        extractCodingBases(refGenome, reqAminoAcids);

        EpitopeUtils.adjustCodingBasesForStrand(this);

        boolean isPhased = phaseMatched();

        int novelUpstreamBases = NovelBaseIndex[FS_UP];
        int novelDownstreamBases = NovelBaseIndex[FS_DOWN];

        if(novelUpstreamBases < 0 || novelDownstreamBases < 0)
        {
            NE_LOGGER.error("ne({}) invalid upBases({} noveBases={}) or downBases({} noveBases={})",
                    this, CodingBases[FS_UP], novelUpstreamBases, CodingBases[FS_DOWN], novelDownstreamBases);
            return;
        }

        // if upstream ends on a phase other than 0, need to take the bases from the downstream gene to make a novel codon
        if(novelUpstreamBases > 0 || novelDownstreamBases > 0 || !isPhased)
        {
            String upstreamBases = CodingBases[FS_UP];
            String downstreamBases = CodingBases[FS_DOWN];

            if(novelUpstreamBases > upstreamBases.length() || novelDownstreamBases > downstreamBases.length())
            {
                NE_LOGGER.debug("ne({}) invalid upBases({} noveBases={}) or downBases({} noveBases={})",
                        this, upstreamBases, novelUpstreamBases, downstreamBases, novelDownstreamBases);
                Valid = false;
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

            CodingBases[FS_UP] = upstreamBases;
            CodingBases[FS_DOWN] = downstreamBases;
        }

        // check superflous upstream bases from longer inserts or INDELs
        int requiredLength = reqAminoAcids * 3;

        if(CodingBases[FS_UP].length() > requiredLength)
            CodingBases[FS_UP] = CodingBases[FS_UP].substring(CodingBases[FS_UP].length() - requiredLength);

        NE_LOGGER.trace("ne({}) upBases({}) novelCodon({}) downBases({}) downNmdBases({})",
                this, CodingBases[FS_UP], EpitopeUtils.checkTrimBases(NovelCodonBases),
                EpitopeUtils.checkTrimBases(CodingBases[FS_DOWN]), NmdBasesMin);
    }

    public void setNonsenseMediatedDecay()
    {
        NmdBasesMin = NmdBasesMax = EpitopeUtils.calcNonMediatedDecayBases(this);

        // also record the bases from the start codon to the stop codon (or end of transcript if there is none)
        int upstreamCodingBases = EpitopeUtils.calcStartCodonBases(this);
        int downstreamCodingBases = EpitopeUtils.calcStopCodonBases(this);

        CodingBasesLengthMin = CodingBasesLengthMax = upstreamCodingBases + downstreamCodingBases;
    }

    public static String trimIncompleteCodons(final String codingBases)
    {
        // strip leading incomplete codon bases from start of coding bases
        int remainder = codingBases.length() % 3;
        if(remainder == 0)
            return codingBases;

        return codingBases.substring(remainder);
    }

    public void setAminoAcids(final RefGenomeInterface refGenome, int reqWildtypeAminoAcids)
    {
        UpstreamAcids = EpitopeUtils.getAminoAcids(trimIncompleteCodons(CodingBases[FS_UP]), false);
        NovelAcid = EpitopeUtils.getAminoAcids(NovelCodonBases, true);
        DownstreamAcids = EpitopeUtils.getAminoAcids(CodingBases[FS_DOWN], true);

        // truncate the novel bases and AAs if a stop codon is encountered
        if(NovelAcid.endsWith(STOP_SYMBOL))
        {
            int novelBases = NovelAcid.length() * 3;
            NovelCodonBases = NovelCodonBases.substring(0, novelBases);
            DownstreamAcids = "";
        }

        // cache bases downstream of mutation in the reference to check for non-novel AAs later on
        int upstreamAAPosStart = orientation(FS_UP) == POS_ORIENT ? ExtPositions[FS_UP][SE_START] : ExtPositions[FS_UP][SE_END];
        byte wtOrient = orientation(FS_UP) == POS_ORIENT ? NEG_ORIENT : POS_ORIENT;
        int requiredBases = CodingBases[FS_UP].length() + reqWildtypeAminoAcids * 3;

        CodingBaseExcerpt wildTypeUpExcerpt = EpitopeUtils.getDownstreamCodingBaseExcerpt(
                refGenome, TransData[FS_UP], chromosome(FS_UP), upstreamAAPosStart, wtOrient,
                requiredBases, true, false, false);

        if(wildTypeUpExcerpt != null)
        {
            String upWildtypeBases = wildTypeUpExcerpt.Bases;

            if(strand(FS_UP) == NEG_STRAND)
                upWildtypeBases = reverseComplementBases(upWildtypeBases);

            UpstreamWildTypeAcids = EpitopeUtils.getAminoAcids(trimIncompleteCodons(upWildtypeBases), true);
        }

        NE_LOGGER.trace("ne({}) upAA({}) novel({}) downAA({})",
                this, UpstreamAcids, EpitopeUtils.checkTrimBases(NovelAcid), EpitopeUtils.checkTrimBases(DownstreamAcids));

        checkStopLost(refGenome, reqWildtypeAminoAcids);
    }

    public abstract void checkStopLost(final RefGenomeInterface refGenome, int reqWildtypeAminoAcids);

    public abstract void setSkippedSpliceSites(final EnsemblDataCache geneTransCache);

    public String aminoAcidString() { return UpstreamAcids + NovelAcid + DownstreamAcids; }

    public NeoEpitopeFile toFile(final int neId, final Set<String> upTransNames, final Set<String> downTransNames)
    {
        return new NeoEpitopeFile(
                neId, variantType(), variantInfo(), variantCopyNumber(), copyNumber(), subclonalLikelihood(),
                TransData[FS_UP].GeneId, TransData[FS_DOWN].GeneId, geneName(FS_UP), geneName(FS_DOWN),
                chromosome(FS_UP), chromosome(FS_DOWN), orientation(FS_UP), orientation(FS_DOWN),
                UpstreamAcids, DownstreamAcids, NovelAcid, NmdBasesMin, NmdBasesMax, CodingBasesLengthMin, CodingBasesLengthMax,
                unsplicedDistance(), skippedDonors(), skippedAcceptors(),
                transcriptsToStr(upTransNames), transcriptsToStr(downTransNames), wildtypeAcids(),
                ExtPositions[FS_UP][SE_START], ExtPositions[FS_UP][SE_END], ExtCodingBases[FS_UP], ExtCigars[FS_UP].toString(),
                ExtPositions[FS_DOWN][SE_START], ExtPositions[FS_DOWN][SE_END], ExtCodingBases[FS_DOWN],
                ExtCigars[FS_DOWN] != null ? ExtCigars[FS_DOWN].toString() : "");
    }
}
