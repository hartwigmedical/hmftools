package com.hartwig.hmftools.compar.resource;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_DIFFERENCES;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_MISMATCH_TYPE;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_NEW_VALUES;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_OLD_VALUES;
import static com.hartwig.hmftools.compar.common.CommonUtils.formOutputFilePrefix;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.MismatchType.DUPLICATE_KEY;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

public class ResourceMismatchWriter implements AutoCloseable
{
    private final BufferedWriter mWriter;

    private static final String FLD_FILE = "File";

    public ResourceMismatchWriter(final ResourceComparConfig config) throws IOException
    {
        mWriter = createBufferedWriter(formOutputFilePrefix(config.OutputDir, "compar_resources", config.OutputId) + TSV_EXTENSION);
        writeRow(FLD_CATEGORY, FLD_FILE, FLD_MISMATCH_TYPE, FLD_ITEM_KEY, FLD_DIFFERENCES, FLD_OLD_VALUES, FLD_NEW_VALUES);
    }

    public void writeMismatches(
            final ResourceFilePair filePair, final List<ResourceMismatch> mismatches, @Nullable final ResourceFile oldFile,
            @Nullable final ResourceFile newFile) throws IOException
    {
        String category = filePair.definition().Type.toString();

        for(ResourceMismatch mismatch : mismatches)
        {
            String type = mismatch.type().toString();
            String diffs = String.join(ITEM_DELIM, mismatch.diffValues());

            if(mismatch.type() == DUPLICATE_KEY)
            {
                for(String oldLine : mismatch.oldLines())
                {
                    writeRow(category, filePair.fileName(), type, mismatch.key(), diffs, rowValues(oldFile, oldLine), "");
                }

                for(String newLine : mismatch.newLines())
                {
                    writeRow(category, filePair.fileName(), type, mismatch.key(), diffs, "", rowValues(newFile, newLine));
                }
            }
            else
            {
                String oldValues = mismatch.oldLines().isEmpty() ? "" : rowValues(oldFile, mismatch.oldLines().get(0));
                String newValues = mismatch.newLines().isEmpty() ? "" : rowValues(newFile, mismatch.newLines().get(0));
                writeRow(category, filePair.fileName(), type, mismatch.key(), diffs, oldValues, newValues);
            }
        }
    }

    @Override
    public void close() throws IOException
    {
        mWriter.close();
    }

    private void writeRow(final String... values) throws IOException
    {
        mWriter.write(String.join(TSV_DELIM, values));
        mWriter.newLine();
    }

    private static String rowValues(final ResourceFile file, final String line)
    {
        return file.fieldValues(line).values().stream()
                .map(x -> format("%s=%s", x.name(), x.displayValue()))
                .collect(Collectors.joining(ITEM_DELIM));
    }
}
