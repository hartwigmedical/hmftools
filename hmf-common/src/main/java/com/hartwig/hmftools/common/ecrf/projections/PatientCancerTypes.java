package com.hartwig.hmftools.common.ecrf.projections;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.google.common.collect.Lists;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientCancerTypes {
    private enum Header {
        cpctId,
        cancerType,
        cancerSubtype
    }

    @NotNull
    abstract String cpctId();

    @NotNull
    abstract String cancerType();

    @NotNull
    abstract String cancerSubtype();

    @NotNull
    private List<String> csvRecord() {
        return Lists.newArrayList(cpctId(), cancerType(), cancerSubtype());
    }

    public static void writeRecords(@NotNull final String outputPath, @NotNull final List<PatientCancerTypes> patientCancerTypes)
            throws IOException {
        final CSVFormat format = CSVFormat.DEFAULT.withNullString("").withHeader(Header.class);
        final CSVPrinter printer = new CSVPrinter(new FileWriter(outputPath), format);
        printer.printRecords(patientCancerTypes.stream().map(PatientCancerTypes::csvRecord).collect(Collectors.toList()));
        printer.close();
    }

    public static List<PatientCancerTypes> readRecords(@NotNull final String filePath) throws IOException {
        final CSVParser parser = CSVParser.parse(new File(filePath),
                Charset.defaultCharset(),
                CSVFormat.DEFAULT.withHeader(Header.class).withSkipHeaderRecord());
        return StreamSupport.stream(parser.spliterator(), false)
                .map(record -> ImmutablePatientCancerTypes.of(record.get(Header.cpctId),
                        record.get(Header.cancerType),
                        record.get(Header.cancerSubtype)))
                .collect(Collectors.toList());
    }
}
