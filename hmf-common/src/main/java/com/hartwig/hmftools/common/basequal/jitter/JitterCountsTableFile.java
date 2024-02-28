package com.hartwig.hmftools.common.basequal.jitter;

import java.io.File;
import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;

public class JitterCountsTableFile
{
    private static String UNIT = "unit";
	private static String NUM_UNITS = "numUnits";

	private static String READ_COUNT = "readCount";
	private static String COUNT_m10 = "count-10";
	private static String COUNT_m9 = "count-9";
	private static String COUNT_m8 = "count-8";
	private static String COUNT_m7 = "count-7";
	private static String COUNT_m6 = "count-6";
	private static String COUNT_m5 = "count-5";
	private static String COUNT_m4 = "count-4";
	private static String COUNT_m3 = "count-3";
	private static String COUNT_m2 = "count-2";
	private static String COUNT_m1 = "count-1";
	private static String COUNT_p0 = "count+0";
	private static String COUNT_p1 = "count+1";
	private static String COUNT_p2 = "count+2";
	private static String COUNT_p3 = "count+3";
	private static String COUNT_p4 = "count+4";
	private static String COUNT_p5 = "count+5";
	private static String COUNT_p6 = "count+6";
	private static String COUNT_p7 = "count+7";
	private static String COUNT_p8 = "count+8";
	private static String COUNT_p9 = "count+9";
	private static String COUNT_p10 = "count+10";

	private static String FILE_EXTENSION = ".ms_table.tsv.gz";

	private static int MAX_UNITS = 20;

	public static String generateFilename(String basePath, String sample)
	{
		return basePath + File.separator + sample + FILE_EXTENSION;
	}

	public static void write(final String filename, @NotNull final Collection<JitterCountsTable> msStatsTables)
	{
		List<String> columns = List.of(UNIT, NUM_UNITS, READ_COUNT,
				COUNT_m10, COUNT_m9, COUNT_m8, COUNT_m7, COUNT_m6, COUNT_m5, COUNT_m4, COUNT_m3, COUNT_m2, COUNT_m1,
				COUNT_p0, COUNT_p1, COUNT_p2, COUNT_p3, COUNT_p4, COUNT_p5, COUNT_p6, COUNT_p7, COUNT_p8, COUNT_p9, COUNT_p10);

		try(DelimFileWriter<JitterCountsTable.Row> writer = new DelimFileWriter<>(filename, columns, (msStatsTableRow, row) ->
			{
				row.set(UNIT,  msStatsTableRow.getRepeatUnit());
				row.set(NUM_UNITS,  msStatsTableRow.refNumUnits);
				row.set(READ_COUNT,  msStatsTableRow.totalReadCount);
				row.set(COUNT_p0, msStatsTableRow.getJitterReadCount(0));

				row.set(COUNT_p10, msStatsTableRow.getJitterReadCount(10));
				row.set(COUNT_p9, msStatsTableRow.getJitterReadCount(9));
				row.set(COUNT_p8, msStatsTableRow.getJitterReadCount(8));
				row.set(COUNT_p7, msStatsTableRow.getJitterReadCount(7));
				row.set(COUNT_p6, msStatsTableRow.getJitterReadCount(6));
				row.set(COUNT_p5, msStatsTableRow.getJitterReadCount(5));
				row.set(COUNT_p4, msStatsTableRow.getJitterReadCount(4));
				row.set(COUNT_p3, msStatsTableRow.getJitterReadCount(3));
				row.set(COUNT_p2, msStatsTableRow.getJitterReadCount(2));
				row.set(COUNT_p1, msStatsTableRow.getJitterReadCount(1));

				row.set(COUNT_m10, msStatsTableRow.getJitterReadCount(-10));
				row.set(COUNT_m9, msStatsTableRow.getJitterReadCount(-9));
				row.set(COUNT_m8, msStatsTableRow.getJitterReadCount(-8));
				row.set(COUNT_m7, msStatsTableRow.getJitterReadCount(-7));
				row.set(COUNT_m6, msStatsTableRow.getJitterReadCount(-6));
				row.set(COUNT_m5, msStatsTableRow.getJitterReadCount(-5));
				row.set(COUNT_m4, msStatsTableRow.getJitterReadCount(-4));
				row.set(COUNT_m3, msStatsTableRow.getJitterReadCount(-3));
				row.set(COUNT_m2, msStatsTableRow.getJitterReadCount(-2));
				row.set(COUNT_m1, msStatsTableRow.getJitterReadCount(-1));
			}))
		{
			for(JitterCountsTable table : msStatsTables)
			{
				for(JitterCountsTable.Row row : table.getRows())
				{
					if(row.refNumUnits <= MAX_UNITS)
					{
						writer.writeRow(row);
					}
				}
			}
		}
	}
}
