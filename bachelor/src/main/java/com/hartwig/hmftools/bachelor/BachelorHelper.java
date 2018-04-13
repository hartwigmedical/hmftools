package com.hartwig.hmftools.bachelor;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.xml.sax.SAXException;

class BachelorHelper {

    private static final Logger LOGGER = LogManager.getLogger(BachelorHelper.class);

    private BachelorHelper() {
    }

    @NotNull
    public static Map<String, Program> loadXML(final Path path) throws IOException, SAXException {
        final BachelorSchema schema = BachelorSchema.make();

        final List<Program> programs = Files.walk(path)
                .filter(p -> p.toString().endsWith(".xml"))
                .map(schema::processXML)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        final Map<String, Program> result = Maps.newHashMap();
        for (final Program p : programs) {
            if (result.containsKey(p.getName())) {
                LOGGER.error("duplicate programs detected: {}", p.getName());
                System.exit(1);
            } else {
                result.put(p.getName(), p);
            }
        }

        return result;
    }
}
