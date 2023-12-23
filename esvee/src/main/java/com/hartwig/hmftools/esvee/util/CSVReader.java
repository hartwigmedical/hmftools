package com.hartwig.hmftools.esvee.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.jetbrains.annotations.Nullable;

public class CSVReader<T> {
    private static final Map<Character, Character> CHARACTER_ESCAPES = Map.of(
            'n', '\n',
            't', '\t',
            'r', '\r',
            '\'', '\'',
            '"', '"',
            '\\', '\\'
    );
    private final Class<? extends T> mItem;
    private final String mFilename;
    private final BufferedReader mReader;
    private final char mDelimiter;
    private int mLineNumber = 1;
    private final List<Setter<T>> mSetters;

    public CSVReader(final Class<? extends T> item, final String filename) throws IOException
    {
        this(item, new File(filename));
    }

    public CSVReader(final Class<? extends T> item, final File file) throws IOException
    {
        mItem = item;
        mFilename = file.getAbsolutePath();
        mReader = new BufferedReader(new FileReader(file));

        @Nullable
        final String header = readLine(); // Skip header
        if (header == null)
        {
            mSetters = List.of();
            mDelimiter = ',';
            return;
        }

        String delimiter = ",";
        String[] parts = header.split(delimiter);
        if (parts.length == 1)
        {
            delimiter = "\t";
            parts = header.split(delimiter);
        }
        mDelimiter = delimiter.charAt(0);

        mSetters = new ArrayList<>();
        for(final String part : parts)
        {
            @Nullable
            final Field field = findField(item, part);
            if (field == null)
            {
                mSetters.add((instance, value) -> {});
                continue;
            }

            field.setAccessible(true);
            if(field.getType() == String.class)
                mSetters.add(field::set);
            else if(field.getType() == Integer.class || field.getType() == int.class)
                mSetters.add((instance, value) -> field.setInt(instance, Integer.parseInt(value)));
            else if(field.getType() == Long.class || field.getType() == long.class)
                mSetters.add((instance, value) -> field.setLong(instance, Long.parseLong(value)));
            else if(field.getType() == Double.class || field.getType() == double.class)
                mSetters.add((instance, value) -> field.setDouble(instance, Double.parseDouble(value)));
            else if(field.getType() == Float.class || field.getType() == float.class)
                mSetters.add((instance, value) -> field.setFloat(instance, Float.parseFloat(value)));
            else if(field.getType() == Boolean.class || field.getType() == boolean.class)
                mSetters.add((instance, value) -> field.setBoolean(instance, StringUtils.parseBoolean(value)));
            else if (field.getType().isEnum())
            {
                final Optional<Method> parseMethod = Arrays.stream(field.getType().getDeclaredMethods())
                        .filter(method -> Modifier.isStatic(method.getModifiers()))
                        .filter(method -> method.getName().equals("parse"))
                        .filter(method -> method.getParameterCount() == 1 && method.getParameterTypes()[0].equals(String.class))
                        .filter(method -> method.getReturnType().equals(field.getType()))
                        .findFirst();
                if (parseMethod.isPresent())
                    mSetters.add((instance, value) -> field.set(instance, parseMethod.get().invoke(instance, value)));
                else
                {
                    //noinspection unchecked,rawtypes
                    mSetters.add((instance, value) -> field.set(instance, Enum.valueOf((Class<Enum>) field.getType(), value)));
                }
            }
            else
                throw new IllegalStateException("");
        }
    }

    @Nullable
    private String readLine() throws IOException
    {
        @Nullable
        final String line = mReader.readLine();
        if (line == null)
            return null;

        return line.replace('\u1680', ' ').replace('\uFEFF', ' ').replace('\uFFFE', ' ').trim();
    }

    @Nullable
    private static Field findField(final Class<?> clazz, final String fieldName)
    {
        try
        {
            return clazz.getDeclaredField(fieldName);
        }
        catch(final NoSuchFieldException ignored)
        {
        }

        return Arrays.stream(clazz.getDeclaredFields())
                .filter(f -> f.getName().equalsIgnoreCase(fieldName))
                .findFirst().orElse(null);
    }

    private T newInstance()
    {
        try
        {
            //noinspection unchecked
            return (T) UnsafeGetter.UNSAFE.allocateInstance(mItem);
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    @Nullable
    public T next() throws IOException
    {
        @Nullable
        final String line = readLine();
        if(line == null)
            return null;

        final List<String> tokens = tokenise(line, mLineNumber++, mDelimiter);
        try
        {
            final T row = newInstance();
            for(int i = 0; i < tokens.size(); i++)
                if (!tokens.get(i).isEmpty())
                {
                    try
                    {
                        mSetters.get(i).setValue(row, tokens.get(i));
                    }
                    catch(final Exception e)
                    {
                        throw new RuntimeException(String.format("Failure while parsing value '%s' on line %s of %s",
                                tokens.get(i), mLineNumber - 1, mFilename), e);
                    }
                }

            return row;
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    public List<T> readToEnd() throws IOException
    {
        final List<T> results = new ArrayList<>();
        while(true)
        {
            final T item = next();
            if(item == null)
                return results;

            results.add(item);
        }
    }

    private static boolean checkCharacter(final String line, final int checkIndex, final char expected) {
        if (line.length() <= checkIndex) {
            return false;
        }
        return line.charAt(checkIndex) == expected;
    }


    private static int takePlainToken(final String line, final int startPosition, final char delimiter, final StringBuilder output) {
        for (int position = startPosition; position < line.length(); position++) {
            final char c = line.charAt(position);
            if (c == delimiter) {
                return position + 1;
            } else {
                output.append(c);
            }
        }
        return line.length();
    }

    private static int advanceToNextToken(final String line, final int lineNumber, final int startPosition, final char delimiter) {
        for (int position = startPosition; position < line.length(); position++) {
            final char c = line.charAt(position);
            if (c == delimiter) {
                return position + 1;
            } else if (!Character.isWhitespace(c)) {
                throw new IllegalStateException("Unexpected character " + c + " while looking for delimiter on line " + lineNumber);
            }
        }
        return line.length();
    }

    /** It is assumed the quote character being used is at `line[startPosition]` */
    private static int takeQuotedToken(final String line, final int lineNumber, final int startPosition, final char delimiter, final StringBuilder output) {
        final char quoteCharacter = line.charAt(startPosition);
        for (int position = startPosition + 1; position < line.length(); position++) {
            final char c = line.charAt(position);

            if (c == quoteCharacter) {
                // Double quote (ex "" or '') inside a string delimited using that character is a valid escape.
                if (checkCharacter(line, position + 1, quoteCharacter)) {
                    output.append(c);
                    position++;
                } else {
                    return advanceToNextToken(line, lineNumber, position, delimiter);
                }
            } else if (c == '\\') {
                if (position + 1 >= line.length()) {
                    throw new IllegalStateException("Unexpected end-of-line while processing escape sequence on line " + lineNumber);
                }
                @Nullable final Character unescaped = CHARACTER_ESCAPES.get(line.charAt(position + 1));
                if (unescaped == null) {
                    throw new IllegalStateException("Illegal escape sequence in string on line "
                            + lineNumber + ": \\" + line.charAt(position + 1));
                }
                output.append(unescaped);
                position++;
            } else {
                output.append(c);
            }
        }

        throw new IllegalStateException("Unclosed string literal beginning at position " + startPosition + " on line " + lineNumber);
    }

    private static List<String> tokenise(final String line, final int lineNumber, final char delimiter) {
        final List<String> tokens = new ArrayList<>();

        final StringBuilder currentToken = new StringBuilder();
        for (int i = 0; i < line.length(); )
        {
            final char c = line.charAt(i++);
            if(Character.isWhitespace(c))
                continue;

            if(c == delimiter)
                tokens.add("");
            else if(c == '"')
            {
                i = takeQuotedToken(line, lineNumber, i - 1, delimiter, currentToken);
                tokens.add(currentToken.toString());
                currentToken.setLength(0);
            }
            else
            {
                i = takePlainToken(line, i - 1, delimiter, currentToken);
                tokens.add(currentToken.toString());
                currentToken.setLength(0);
            }
        }

        return tokens;
    }

    private interface Setter<T> {
        void setValue(T instance, String rawValue) throws Exception;
    }
}
