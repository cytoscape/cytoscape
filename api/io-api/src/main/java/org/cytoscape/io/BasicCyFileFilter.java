package org.cytoscape.io;

/*
 * #%L
 * Cytoscape IO API (io-api)
 * $Id:$
 * $HeadURL:$
 * %%
 * Copyright (C) 2006 - 2013 The Cytoscape Consortium
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by the Free Software Foundation, either version 2.1 of the 
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public 
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-2.1.html>.
 * #L%
 */


import org.cytoscape.io.util.StreamUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URI;
import java.util.HashSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This is meant to be an basic implementation of {@link CyFileFilter} that can
 * either be used directly or extended to provide different acceptance criteria.
 * Only the accepts() methods may be overridden.
 *
 * @CyAPI.Abstract.Class
 * @CyAPI.InModule io-api
 */
public class BasicCyFileFilter implements CyFileFilter {

	protected final Set<String> extensions;
	protected final Set<String> contentTypes;
	protected final String description;
	protected final StreamUtil streamUtil;
	/** Type of data that this filter applies to. */
	protected final DataCategory category;
	private static final Logger logger = LoggerFactory.getLogger(BasicCyFileFilter.class);

	/**
	 * Creates a file filter from the specified arguments. Note that a "."
	 * before the extension is not needed and will be ignored.
	 *
	 * @param extensions
	 *            The set of valid extensions for this filter.
	 * @param contentTypes
	 *            The set of valid MIME content types that this filter should
	 *            recognize.
	 * @param description
	 *            A human readable description of the filter.
	 * @param category
	 *            The type of data this filter is meant to support.
	 * @param streamUtil
	 *            An instance of the StreamUtil service.
	 */
	public BasicCyFileFilter(final Set<String> extensions,
			final Set<String> contentTypes, final String description,
			final DataCategory category, StreamUtil streamUtil) {

		this.extensions = cleanStringSet(extensions);
		this.contentTypes = contentTypes;
		this.category = category;

		final StringBuilder builder = new StringBuilder();
		builder.append(description == null ? "(" : description + " (");

		for (String ex : this.extensions) {
			builder.append("*." + ex + ", ");
		}

		String d = builder.toString();
		d = d.substring(0, d.length() - 2);
		d += ")";

		this.description = d;
		this.streamUtil = streamUtil;
	}

	/**
	 * Creates a file filter from the specified arguments. Note that a "."
	 * before the extension is not needed and will be ignored.
	 *
	 * @param extensions
	 *            The set of valid extensions for this filter.
	 * @param contentTypes
	 *            The set of valid MIME content types that this filter should
	 *            recognize.
	 * @param description
	 *            A human readable description of the filter.
	 * @param category
	 *            The type of data this filter is meant to support.
	 * @param streamUtil
	 *            An instance of the StreamUtil service.
	 */
	public BasicCyFileFilter(final String[] extensions, final String[] contentTypes, final String description,
			final DataCategory category, StreamUtil streamUtil) {
		this(createSet(extensions), createSet(contentTypes), description, category, streamUtil);
	}

	/**
	 * This removes null and empty (white space only) Strings from Set strings
	 *
	 * @param strings
	 *            a Set of Strings to be cleaned
	 * @return a SorteSet of Strings containing only non-empty Strings
	 */
	private final static SortedSet<String> cleanStringSet(final Set<String> strings) {
		final SortedSet<String> cleaned_strings = new TreeSet<String>();
		if (strings != null) {
			for (final String string : strings) {
				if (string != null && string.trim().length() > 0) {
					cleaned_strings.add(string.trim());
				}
			}
		}
		return cleaned_strings;
	}
	
	private static Set<String> createSet(String[] values) {
		Set<String> set = new HashSet<String>();
		for (String v : values)
			set.add(v);
		return set;
	}

	/**
	 * {@inheritDoc}
	 */
	public boolean accepts(URI uri, DataCategory category) {

		// Check data category
		if (category != this.category)
			return false;

		if (extensionsMatch(uri))
			return true;
		else
			return false;

	}

	private boolean extensionsMatch(URI uri) {
		final String extension = getExtension(uri.toString());
		//extension is never null anymore, but can be empty string, which works ok
		if (extensions.contains(extension))
			return true;
		else
			return false;
	}

	/**
	 * This method always returns false in this particular implementation. You
	 * must extend this class and override this method to get alternative
	 * behavior. Ideally this method would return true if this class is capable
	 * of processing the specified InputStream.
	 *
	 * @param stream
	 *            The stream that references the file we'd like to read.
	 * @param category
	 *            The type of input that we're considering.
	 * @return Always returns false in this particular implementation.
	 */
	public boolean accepts(InputStream stream, DataCategory category) {
		return false;
	}

	/**
	 * {@inheritDoc}
	 */
	public final Set<String> getExtensions() {
		return extensions;
	}

	/**
	 * {@inheritDoc}
	 */
	public final Set<String> getContentTypes() {
		return contentTypes;
	}

	/**
	 * {@inheritDoc}
	 */
	public final String getDescription() {
		return description;
	}

	/**
	 * {@inheritDoc}
	 */
	public final DataCategory getDataCategory() {
		return category;
	}

	/**
	 * Returns a human readable description of this class.
	 *
	 * @return a human readable description of this class.
	 */
	@Override
	public String toString() {
		final StringBuilder builder = new StringBuilder();
		builder.append(description + " [category: " + category + "]  [extensions: ");

		for (String ext : extensions)
			builder.append(ext + ",");

		String s = builder.toString();
		s += "]   [contentTypes: ";
		for (String c : contentTypes)
			s += c + ",";
		s += "]";

		return s;
	}

	/**
	 * Returns a string of the characters following the last '.' in the input
	 * string, which is to say the file extension assuming that the input string
	 * represents a file name. Will return null if no '.' is found.
	 *
	 * @param filename
	 *            the file name as a string.
	 * @return a string representing the file extension of the input string.
	 */
	protected final String getExtension(String filename) {
		return CommonsIOFilenameUtils.getExtension(filename).toLowerCase();
	}

	/**
	 * Returns a string containing the specified number of lines from the
	 * beginning of the file. This is useful for testing input streams.
	 *
	 * @param stream
	 *            the input stream from which to read the header.
	 * @param numLines
	 *            the number of lines from the beginning of the file.
	 * @return a string containing the specified number of lines from the
	 *         beginning of the file.
	 */
	protected final String getHeader(InputStream stream, int numLines) {

		String header;
		BufferedReader br = new BufferedReader(new InputStreamReader(stream));

		try {
			header = parseHeader(br, numLines);
		} catch (IOException ioe) {
			logger.warn("failed to read header from stream", ioe);
			header = "";
		} finally {
			if (br != null)
				try {
					br.close();
				} catch (IOException e) {
				}

			br = null;
		}

		return header;
	}

	private final String parseHeader(BufferedReader bufferedReader, int numLines) throws IOException {
		StringBuilder header = new StringBuilder();

		try {
			String line = bufferedReader.readLine();

			int lineCount = 0;

			while ((line != null) && (lineCount < numLines)) {
				header.append(line + "\n");
				line = bufferedReader.readLine();
				lineCount++;
			}
		} finally {
			if (bufferedReader != null)
				bufferedReader.close();
		}

		return header.toString();
	}

	/*
	 * Licensed to the Apache Software Foundation (ASF) under one or more
	 * contributor license agreements.  See the NOTICE file distributed with
	 * this work for additional information regarding copyright ownership.
	 * The ASF licenses this file to You under the Apache License, Version 2.0
	 * (the "License"); you may not use this file except in compliance with
	 * the License.  You may obtain a copy of the License at
	 * 
	 *      http://www.apache.org/licenses/LICENSE-2.0
	 * 
	 * Unless required by applicable law or agreed to in writing, software
	 * distributed under the License is distributed on an "AS IS" BASIS,
	 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	 * See the License for the specific language governing permissions and
	 * limitations under the License.
	 */
	/**
	 * Taken from org.apache.commons.io.FilenameUtils 2.1 to avoid
	 * adding dependency on commons-io.
	 */
	private static class CommonsIOFilenameUtils {
	    /**
	     * The extension separator character.
	     * @since Commons IO 1.4
	     */
	    public static final char EXTENSION_SEPARATOR = '.';

	    /**
	     * The Unix separator character.
	     */
	    private static final char UNIX_SEPARATOR = '/';

	    /**
	     * The Windows separator character.
	     */
	    private static final char WINDOWS_SEPARATOR = '\\';
	    
	    /**
	     * Gets the extension of a filename.
	     * <p>
	     * This method returns the textual part of the filename after the last dot.
	     * There must be no directory separator after the dot.
	     * <pre>
	     * foo.txt      --> "txt"
	     * a/b/c.jpg    --> "jpg"
	     * a/b.txt/c    --> ""
	     * a/b/c        --> ""
	     * </pre>
	     * <p>
	     * The output will be the same irrespective of the machine that the code is running on.
	     *
	     * @param filename the filename to retrieve the extension of.
	     * @return the extension of the file or an empty string if none exists or <code>null</code>
	     * if the filename is <code>null</code>.
	     */
	    public static String getExtension(String filename) {
	        if (filename == null) {
	            return null;
	        }
	        int index = indexOfExtension(filename);
	        if (index == -1) {
	            return "";
	        } else {
	            return filename.substring(index + 1);
	        }
	    }

	    /**
	     * Returns the index of the last extension separator character, which is a dot.
	     * <p>
	     * This method also checks that there is no directory separator after the last dot.
	     * To do this it uses {@link #indexOfLastSeparator(String)} which will
	     * handle a file in either Unix or Windows format.
	     * <p>
	     * The output will be the same irrespective of the machine that the code is running on.
	     * 
	     * @param filename  the filename to find the last path separator in, null returns -1
	     * @return the index of the last separator character, or -1 if there
	     * is no such character
	     */
	    public static int indexOfExtension(String filename) {
	        if (filename == null) {
	            return -1;
	        }
	        int extensionPos = filename.lastIndexOf(EXTENSION_SEPARATOR);
	        int lastSeparator = indexOfLastSeparator(filename);
	        return (lastSeparator > extensionPos ? -1 : extensionPos);
	    }
	    
	    /**
	     * Returns the index of the last directory separator character.
	     * <p>
	     * This method will handle a file in either Unix or Windows format.
	     * The position of the last forward or backslash is returned.
	     * <p>
	     * The output will be the same irrespective of the machine that the code is running on.
	     * 
	     * @param filename  the filename to find the last path separator in, null returns -1
	     * @return the index of the last separator character, or -1 if there
	     * is no such character
	     */
	    public static int indexOfLastSeparator(String filename) {
	        if (filename == null) {
	            return -1;
	        }
	        int lastUnixPos = filename.lastIndexOf(UNIX_SEPARATOR);
	        int lastWindowsPos = filename.lastIndexOf(WINDOWS_SEPARATOR);
	        return Math.max(lastUnixPos, lastWindowsPos);
	    }
	}
}
