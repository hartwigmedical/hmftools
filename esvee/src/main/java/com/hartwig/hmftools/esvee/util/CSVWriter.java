package com.hartwig.hmftools.esvee.util;

import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.annotation.Nullable;

import org.apache.commons.lang3.tuple.Pair;

public enum CSVWriter
{
    ;

    public static <T> void writeCSV(String filename, Class<? extends T> type, Stream<T> values)
    {
        writeCSV(filename, type, values::iterator);
    }

    public static <T> void writeCSV(String filename, Class<? extends T> type, Iterable<T> values)
    {
        try
        {
            Files.writeString(Path.of(filename), writeCSV(new StringBuilder(), type, values));
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static <T> String writeCSV(Class<? extends T> type, Iterable<T> values)
    {
        return writeCSV(new StringBuilder(), type, values).toString();
    }

    public static <T> StringBuilder writeCSV(StringBuilder sb, Class<? extends T> type, Iterable<T> values)
    {
        final var context = createContext(type, ',');
        appendCSVHeader(context, sb);
        values.forEach(value -> appendCSVValue(context, sb, value));
        return sb;
    }

    private static <T> void appendCSVHeader(CSVWriterContext context, StringBuilder sb)
    {
        sb.append(context.FieldNames.stream()
                .collect(Collectors.joining(String.valueOf(context.FieldDelimiter))));
        sb.append("\n");
    }

    private static <T> void appendCSVValue(CSVWriterContext context, StringBuilder sb, T value)
    {
        sb.append(context.FieldGetters.stream()
                .map(getter -> valueToString(context, getter.apply(value)))
                .collect(Collectors.joining(String.valueOf(context.FieldDelimiter))));
        sb.append("\n");
    }

    private static String valueToString(CSVWriterContext context, @Nullable Object value)
    {
        if(value == null)
            return "";
        else if(value instanceof String && needsEscaping(context.FieldDelimiter, (String) value))
            return String.format("\"%s\"", ((String) value).replaceAll("\"", "\"\""));
        else if(value instanceof CSVValue)
            return valueToString(context, ((CSVValue) value).asValueForCSV());
        else if(value instanceof Iterable)
            return StreamSupport.stream(((Iterable<?>) value).spliterator(), false)
                    .map(innerValue -> valueToString(context, innerValue))
                    .collect(Collectors.joining(";"));
        else
            return value.toString();
    }

    private static boolean needsEscaping(char separator, String value)
    {
        return value.startsWith("\"") || value.contains(String.valueOf(separator));
    }

    private static final Set<Class<?>> BASIC_TYPES = Set.of(
            byte.class, short.class, int.class, long.class, float.class, double.class,
            String.class,
            boolean.class, Boolean.class);

    private static CSVWriterContext createContext(Class<?> valueType, char delimiter)
    {
        var pair = createGetters(valueType);
        return new CSVWriterContext(valueType, delimiter, pair.getLeft(), pair.getRight());
    }

    private static Function<Object, Object> getter(Field field) {
        field.setAccessible(true);
        return instance ->
        {
            try
            {
                return field.get(instance);
            }
            catch(IllegalAccessException e)
            {
                throw new RuntimeException(e);
            }
        };
    }

    private static Pair<List<String>, List<Function<Object, Object>>> createGetters(Class<?> type)
    {
        final List<String> fieldNames = new ArrayList<>();
        final List<Function<Object, Object>> fieldGetters = new ArrayList<>();
        for(final Field field : type.getDeclaredFields())
        {
            if(Modifier.isStatic(field.getModifiers()) || field.getName().contains("__"))
                continue;

            if(BASIC_TYPES.contains(field.getType()) || Number.class.isAssignableFrom(field.getType())
                    || field.getType().isEnum() || Iterable.class.isAssignableFrom(field.getType())
                    || CSVValue.class.isAssignableFrom(field.getType()))
            {
                fieldNames.add(field.getName());
                fieldGetters.add(getter(field));
            }
            else
            {
                final Function<Object, Object> rootGetter = getter(field);
                final var pair = createGetters(field.getType());
                for(int i = 0; i < pair.getLeft().size(); i++)
                {
                    fieldNames.add(field.getName() + "_" + pair.getLeft().get(i));
                    fieldGetters.add(rootGetter.andThen(pair.getRight().get(i)));
                }
            }
        }

        return Pair.of(fieldNames, fieldGetters);
    }

    private static class CSVWriterContext
    {
        public final Class<?> ValueType;
        public final char FieldDelimiter;
        public final List<String> FieldNames;
        public final List<Function<Object, Object>> FieldGetters;

        private CSVWriterContext(final Class<?> valueType, final char fieldDelimiter, final List<String> fieldNames,
                final List<Function<Object, Object>> fieldGetters)
        {
            ValueType = valueType;
            FieldDelimiter = fieldDelimiter;
            FieldNames = fieldNames;
            FieldGetters = fieldGetters;
        }
    }

    public interface CSVValue {
        Object asValueForCSV();
    }
}
