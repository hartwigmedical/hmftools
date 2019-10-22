package com.hartwig.hmftools.linx.neoepitope;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class KmerGenerator
{
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final RefGenomeInterface mRefGenome;
    private final String mOutputDir;
    private BufferedWriter mFileWriter;

    private static final Logger LOGGER = LogManager.getLogger(KmerGenerator.class);

    public KmerGenerator(final RefGenomeInterface refGenome, final SvGeneTranscriptCollection geneTransCache, final String outputDir)
    {
        mGeneTransCache = geneTransCache;
        mRefGenome = refGenome;

        mOutputDir = outputDir;
        mFileWriter = null;
    }


    private TranscriptData getTranscriptData(final Transcript transcript)
    {
        final TranscriptData transData = mGeneTransCache.getTranscriptData(transcript.gene().StableId, transcript.StableId);

        if(transData == null)
        {
            LOGGER.error("gene({}) transcript({}) data not found", transcript.gene().GeneName, transcript.StableId);
            return null;
        }

        return transData;
    }

    private String getBaseString(final String chromosome, long posStart, long posEnd)
    {
        String baseString = mRefGenome.getBaseString(chromosome, posStart, posEnd);

        return baseString;
    }

    private void writeData(final String sampleId)
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir + "LNX_KMER_STRINGS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("Chromosome,PosStart,PosEnd,BaseString");

                mFileWriter.newLine();
            }

            /*
            mFileWriter.write(String.format("%s,%s,%s",
                    sampleId, fusion.name(), fusion.upstreamTrans().geneName().equals(fusion.downstreamTrans().geneName())));

            mFileWriter.write(String.format(",%s,%s,%s,%d",
                    data.upstreamAcids(), data.downstreamAcids(), data.novelAcid(), data.downstreamNmdBases()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                final GeneAnnotation gene = trans.gene();

                mFileWriter.write(String.format(",%d,%s,%d,%d",
                        gene.id(), gene.chromosome(), gene.position(), gene.orientation()));

                mFileWriter.write(String.format(",%s,%d,%s,%s,%d,%d",
                        trans.StableId, gene.Strand, trans.regionType(), trans.codingType(),
                        trans.nextSpliceExonRank(), trans.nextSpliceExonPhase()));
            }
            */

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing k-mer output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }



}
