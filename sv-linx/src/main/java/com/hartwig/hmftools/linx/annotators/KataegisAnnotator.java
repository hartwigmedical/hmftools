package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvBreakend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class KataegisAnnotator
{
    private Map<String,Map<String,List<KataegisData>>> mSampleChrData;
    private final String mOutputDir;
    private BufferedWriter mFileWriter;

    private static int PROXIMITY_THRESHOLD = 100000;

    private static final Logger LOGGER = LogManager.getLogger(KataegisAnnotator.class);

    public KataegisAnnotator(final String outputDir)
    {
        mOutputDir = outputDir;
        mSampleChrData = Maps.newHashMap();
        mFileWriter = null;
    }

    public void annotateVariants(final String sampleId, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        final Map<String,List<KataegisData>> sampleData = mSampleChrData.get(sampleId);

        if(sampleData == null || sampleData.isEmpty())
            return;

        for(final Map.Entry<String,List<KataegisData>> entry : sampleData.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<KataegisData> katDataList = entry.getValue();

            if(katDataList.isEmpty())
                continue;

            final List<SvBreakend> breakendList = chrBreakendMap.get(chromosome);

            if(breakendList == null)
            {
                katDataList.forEach(x -> writeSvData(sampleId, x,null));
            }
            else
            {
                for(final KataegisData katData : katDataList)
                {
                    SvBreakend closestBreakend = findClosestBreakend(katData, breakendList);
                    writeSvData(sampleId, katData, closestBreakend);
                }
            }
        }
    }

    private SvBreakend findClosestBreakend(final KataegisData katData, final List<SvBreakend> breakendList)
    {
        // find the nearest facing breakend
        SvBreakend closestBreakend = null;
        long shortestDistance = 0;

        for (final SvBreakend breakend : breakendList)
        {
            if (breakend.position() < katData.PosStart)
            {
                if (breakend.orientation() == 1)
                    continue;

                if (katData.PosStart - breakend.position() > PROXIMITY_THRESHOLD)
                    continue;

                closestBreakend = breakend;
                shortestDistance = katData.PosStart - breakend.position();
            }
            else if (breakend.position() <= katData.PosEnd)
            {
                continue;
            }
            else
            {
                long distance = breakend.position() - katData.PosEnd;

                if (distance > PROXIMITY_THRESHOLD)
                    break;

                if (breakend.orientation() == -1)
                    continue;

                if (closestBreakend == null || shortestDistance > distance)
                {
                    closestBreakend = breakend;
                }

                break;
            }
        }

        return closestBreakend;
    }

    private void writeSvData(final String sampleId, final KataegisData data, final SvBreakend breakend)
    {
        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir + "SVA_KATAEGIS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Chromosome,PosStart,PosEnd,KataegisId,SnvCount");
                mFileWriter.write(",SvId,SvPosition,SvIsStart,SvOrient,Distance");
                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%s,%d,%d,%s,%d",
                    sampleId, data.Chromosome, data.PosStart, data.PosStart, data.Id, data.SnvCount));

            if(breakend != null)
            {
                long distance = min(abs(data.PosStart - breakend.position()), abs(data.PosEnd - breakend.position()));

                mFileWriter.write(String.format(",%d,%d,%s,%d,%d",
                        breakend.getSV().id(), breakend.position(), breakend.usesStart(), breakend.orientation(), distance));
            }
            else
            {
                mFileWriter.write(",-1,,,0,-0");
            }

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing kataegis output file: {}", e.toString());
        }
    }

    private final static int CSV_REQUIRED_FIELDS = 6;

    public void loadKataegisData(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int recordCount = 0;

            // SampleId,KataegisId,SnvCount,Chromosome,MinPosition,MaxPosition

            String line;
            String currentSample = "";
            Map<String,List<KataegisData>> chrData = null;

            while ((line = fileReader.readLine()) != null)
            {
                if(line.contains("SampleId"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < CSV_REQUIRED_FIELDS)
                    continue;

                final String sampleId = items[0];
                final String chromosome = items[3];

                KataegisData data = new KataegisData(
                        chromosome,
                        Long.parseLong(items[4]),
                        Long.parseLong(items[5]),
                        items[1],
                        Integer.parseInt(items[2]));

                if(!currentSample.equals(sampleId))
                {
                    currentSample = sampleId;
                    chrData = Maps.newHashMap();
                    mSampleChrData.put(sampleId, chrData);
                }

                List<KataegisData> dataList = chrData.get(chromosome);

                if(dataList == null)
                {
                    dataList = Lists.newArrayList();
                    chrData.put(chromosome, dataList);
                }

                dataList.add(data);
                ++recordCount;
            }

            LOGGER.debug("loaded {} kataegis data records, samples({})", recordCount, mSampleChrData.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read kataegis CSV file({})", filename);
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

}
