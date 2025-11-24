package com.hartwig.hmftools.esvee.caller;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.fromByteStr;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.output.BreakendWriter.FLD_BREAKEND_INS_SEQ;
import static com.hartwig.hmftools.esvee.assembly.output.BreakendWriter.FLD_BREAKEND_MATE_CHR;
import static com.hartwig.hmftools.esvee.assembly.output.BreakendWriter.FLD_BREAKEND_MATE_ORIENT;
import static com.hartwig.hmftools.esvee.assembly.output.BreakendWriter.FLD_BREAKEND_MATE_POSITION;
import static com.hartwig.hmftools.esvee.assembly.output.BreakendWriter.FLD_SV_TYPE;
import static com.hartwig.hmftools.esvee.common.CommonUtils.compareJunctions;
import static com.hartwig.hmftools.esvee.common.FileCommon.ESVEE_FILE_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.formEsveeInputFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.formOutputFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.WriteType;

public class BreakendWriter
{
    private final CallerConfig mConfig;

    public BreakendWriter(final CallerConfig config)
    {
        mConfig = config;
    }

    public void writeBreakends(final SvDataCache svDataCache)
    {
        String esveeInputDir = pathFromFile(mConfig.VcfFile);

        String initialBreakendTsv = formEsveeInputFilename(
                esveeInputDir, mConfig.fileSampleId(), WriteType.BREAKEND.fileId(), mConfig.OutputId);

        if(!Files.exists(Paths.get(initialBreakendTsv)) && mConfig.OutputId != null)
        {
            initialBreakendTsv = formEsveeInputFilename(
                    esveeInputDir, mConfig.fileSampleId(), WriteType.BREAKEND.fileId(), null);
        }

        if(!Files.exists(Paths.get(initialBreakendTsv)))
        {
            SV_LOGGER.warn("cannot write breakend TSV since input TSV({}) not found", initialBreakendTsv);
            return;
        }

        String newBreakendTsv = formOutputFile(
                mConfig.OutputDir, mConfig.fileSampleId(), ESVEE_FILE_ID, "final_breakend.tsv", mConfig.OutputId);

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(initialBreakendTsv));

            String headerLine = fileReader.readLine();

            Map<String,List<BreakendData>> initialBreakendMap = loadBreakendData(fileReader, headerLine);

            if(initialBreakendMap.isEmpty())
                return;

            SV_LOGGER.debug("annotating {} existing breakends", initialBreakendMap.values().stream().mapToInt(x -> x.size()).sum());

            BufferedWriter writer = createBufferedWriter(newBreakendTsv);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(headerLine);

            // additional fields
            sj.add("Germline");
            sj.add("Filter");
            sj.add("VAF");
            sj.add("PONCount");
            writer.write(sj.toString());
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                List<Breakend> breakendList = svDataCache.getBreakendMap().get(chrStr);
                List<BreakendData> breakendDataList = initialBreakendMap.get(chrStr);

                if(breakendList == null || breakendDataList == null)
                    continue;

                int dataIndex = 0;

                for(Breakend breakend : breakendList)
                {
                    if(dataIndex >= breakendDataList.size())
                        break;

                    BreakendData breakendData = breakendDataList.get(dataIndex);

                    if(!breakendData.matches(breakend))
                    {
                        // search forward and backwards looking for a match
                        boolean matchFound = false;

                        for(int d = 0; d <= 1; ++d)
                        {
                            boolean searchBack = (d == 0);

                            int testDataIndex = searchBack ? dataIndex - 1 : dataIndex + 1;

                            while(testDataIndex >= 0 && testDataIndex < breakendDataList.size())
                            {
                                breakendData = breakendDataList.get(testDataIndex);

                                if(breakendData.Position > breakend.Position && !searchBack)
                                    break;
                                else if(breakendData.Position < breakend.Position && searchBack)
                                    break;

                                if(breakendData.matches(breakend))
                                {
                                    dataIndex = testDataIndex;
                                    matchFound = true;
                                    break;
                                }

                                if(searchBack)
                                    --testDataIndex;
                                else
                                    ++testDataIndex;
                            }

                            if(matchFound)
                                break;
                        }
                    }

                    if(breakendData.matches(breakend))
                    {
                        writeBreakendData(writer, breakend, breakendData);
                        ++dataIndex;
                    }
                    else
                    {
                        SV_LOGGER.debug("unmatched breakend({})", breakend);
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to rewrite breakends: {}", e.toString());
        }
    }

    private void writeBreakendData(final BufferedWriter writer, final Breakend breakend, final BreakendData breakendData) throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(breakendData.ValuesStr);

        sj.add(String.valueOf(breakend.sv().isGermline()));

        Set<FilterType> filters = breakend.sv().filters();
        String filtersStr = filters.isEmpty() ? PASS : filters.stream().map(x -> x.vcfTag()).collect(Collectors.joining(ITEM_DELIM));

        sj.add(filtersStr);

        sj.add(format("%.3f", breakend.calcAllelicFrequency()));

        sj.add(String.valueOf(breakend.sv().ponCount()));

        writer.write(sj.toString());
        writer.newLine();
    }

    private class BreakendData implements Comparable<BreakendData>
    {
        public final String BreakendId;
        public final String Chromosome;
        public final int Position;
        public final Orientation Orient;
        public final StructuralVariantType Type;
        public final String MateChromosome;
        public final int MatePosition;
        public final Orientation MateOrient;
        public final String InsertSequence;

        public final String ValuesStr;

        public BreakendData(
                final String breakendId, final String chromosome, final int position, final Orientation orient, final StructuralVariantType type,
                final String mateChromosome, final int matePosition, final Orientation mateOrient,
                final String insertSequence, final String valuesStr)
        {
            BreakendId = breakendId;
            Chromosome = chromosome;
            Position = position;
            Orient = orient;
            Type = type;
            MateChromosome = mateChromosome;
            MatePosition = matePosition;
            MateOrient = mateOrient;
            InsertSequence = insertSequence;
            ValuesStr = valuesStr;
        }

        public boolean matches(final Breakend other)
        {
            if(!BreakendId.equals(other.VcfId))
                return false;

            if(Type != other.sv().type())
            {
                // allow DUP vs INS differences
                if(!(Type == DUP && other.sv().type() == INS) && !(Type == INS && other.sv().type() == DUP))
                    return false;
            }

            if(!Chromosome.equals(other.Chromosome) || Position != other.Position || Orient != other.Orient)
                return false;

            if(Type == StructuralVariantType.SGL)
                return InsertSequence.equals(other.InsertSequence);

            return MateChromosome.equals(other.otherBreakend().Chromosome)
                    && MatePosition == other.otherBreakend().Position && MateOrient == other.otherBreakend().Orient;
        }

        @Override
        public int compareTo(final BreakendData other)
        {
            return compareJunctions(Chromosome, other.Chromosome, Position, other.Position, Orient, other.Orient);
        }

        public String toString() { return format("%s:%d:%d %s", Chromosome, Position, Orient.asByte(), Type); }
    }

    private Map<String,List<BreakendData>> loadBreakendData(final BufferedReader fileReader, final String headerLine)
    {
        try
        {
            Map<String,List<BreakendData>> chrBreakendDataMap = Maps.newHashMap();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(headerLine, TSV_DELIM);

            int idIndex = fieldsIndexMap.get("Id");
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int orientIndex = fieldsIndexMap.get(FLD_ORIENTATION);
            int mateChrIndex = fieldsIndexMap.get(FLD_BREAKEND_MATE_CHR);
            int matePosIndex = fieldsIndexMap.get(FLD_BREAKEND_MATE_POSITION);
            int mateOrientIndex = fieldsIndexMap.get(FLD_BREAKEND_MATE_ORIENT);
            int typeIndex = fieldsIndexMap.get(FLD_SV_TYPE);
            int insSqeIndex = fieldsIndexMap.get(FLD_BREAKEND_INS_SEQ);

            List<BreakendData> breakendList = null;
            String currentChromosome = "";
            String line;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                String breakendId = values[idIndex];
                String chromosome = values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);
                Orientation orientation = fromByteStr(values[orientIndex]);

                StructuralVariantType svType = StructuralVariantType.valueOf(values[typeIndex]);

                String mateChromosome = "";
                int matePosition = -1;
                Orientation mateOrientation = FORWARD;

                if(svType != StructuralVariantType.SGL)
                {
                    mateChromosome = values[mateChrIndex];
                    matePosition = Integer.parseInt(values[matePosIndex]);
                    mateOrientation = fromByteStr(values[mateOrientIndex]);
                }

                String insSequence = values[insSqeIndex];

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;

                    breakendList = chrBreakendDataMap.get(chromosome);

                    if(breakendList == null)
                    {
                        breakendList = Lists.newArrayList();
                        chrBreakendDataMap.put(chromosome, breakendList);
                    }
                }

                breakendList.add(new BreakendData(
                        breakendId, chromosome, position, orientation, svType, mateChromosome, matePosition, mateOrientation, insSequence, line));
            }

            // sort ahead of comparison
            for(List<BreakendData> breakends : chrBreakendDataMap.values())
            {
                Collections.sort(breakends);
            }

            return chrBreakendDataMap;
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read breakend file: {}", exception.toString());
            return null;
        }
    }
}
