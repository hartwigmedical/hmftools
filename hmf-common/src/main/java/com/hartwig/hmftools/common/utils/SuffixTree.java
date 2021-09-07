package com.hartwig.hmftools.common.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/***
 Taken (and modified) from https://github.com/eugenp/tutorials/tree/master/algorithms-searching under MIT License:

 Copyright (c) 2017 Eugen Paraschiv

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 *
 */

public class SuffixTree
{
    private final Node mRoot;
    private final String mFullText;

    private static final String WORD_TERMINATION = "$";
    private static final int POSITION_UNDEFINED = -1;

    public SuffixTree(String text)
    {
        assert (!text.contains(WORD_TERMINATION));

        mRoot = new Node("", POSITION_UNDEFINED, null);
        for(int i = 0; i < text.length(); i++)
        {
            addSuffix(text.substring(i) + WORD_TERMINATION, i);
        }
        mFullText = text;
    }

    public boolean contains(String pattern)
    {
        return findNode(pattern) != null;
    }

    public int endsWith(String pattern)
    {
        List<Node> nodes = getAllNodesInTraversePath(pattern, mRoot, true);
        if(nodes.isEmpty())
        {
            return 0;
        }

        Node lastNode = nodes.get(nodes.size() - 1);
        if(lastNode == null)
        {
            return 0;
        }

        int parentLength = nodes.stream().limit(nodes.size() - 1).mapToInt(x -> x.text.length()).sum();
        String lastNodeText = lastNode.text.replace(WORD_TERMINATION, "");
        if(parentLength + lastNodeText.length() > pattern.length())
        {
            return 0;
        }

        if(!lastNode.getChildren().isEmpty())
        {
            if(lastNode.getChildren().stream().noneMatch(x -> x.text.equals(WORD_TERMINATION)))
            {
                return 0;
            }
        }

        int remainingLength = Math.min(lastNodeText.length(), pattern.length() - parentLength);
        for(int lastNodeIndex = 0; lastNodeIndex < remainingLength; lastNodeIndex++)
        {
            int patternIndex = lastNodeIndex + parentLength;
            if(pattern.charAt(patternIndex) != lastNodeText.charAt(lastNodeIndex))
            {
                return 0;
            }
        }

        return parentLength + lastNodeText.length();
    }

    public int anyIndexOf(String pattern)
    {
        Node lastNode = findNode(pattern);
        if(lastNode != null)
        {
            return getPositions(lastNode).get(0);
        }

        return -1;
    }

    @NotNull
    public List<Integer> indices(String pattern)
    {
        Node lastNode = findNode(pattern);
        if(lastNode != null)
        {
            return getPositions(lastNode).stream().sorted().collect(Collectors.toList());
        }

        return Collections.emptyList();
    }

    @Nullable
    private Node findNode(String pattern)
    {
        List<Node> nodes = getAllNodesInTraversePath(pattern, mRoot, false);
        if(nodes.size() > 0)
        {
            return nodes.get(nodes.size() - 1);
        }

        return null;
    }

    List<Node> findPartialNode(String pattern)
    {
        return getAllNodesInTraversePath(pattern, mRoot, true);
    }

    @VisibleForTesting
    List<String> markPatternInText(String pattern)
    {
        return indices(pattern).stream().map(m -> markPatternInText(m, pattern)).collect(Collectors.toList());
    }

    private void addSuffix(String suffix, int position)
    {
        List<Node> nodes = getAllNodesInTraversePath(suffix, mRoot, true);
        if(nodes.size() == 0)
        {
            addChildNode(mRoot, suffix, position);
        }
        else
        {
            Node lastNode = nodes.remove(nodes.size() - 1);
            String newText = suffix;
            if(nodes.size() > 0)
            {
                String existingSuffixUptoLastNode = nodes.stream().map(Node::getText).reduce("", String::concat);

                // Remove prefix from newText already included in parent
                newText = newText.substring(existingSuffixUptoLastNode.length());
            }
            extendNode(lastNode, newText, position);
        }
    }

    private List<Integer> getPositions(Node node)
    {
        List<Integer> positions = new ArrayList<>();
        if(node.getText().endsWith(WORD_TERMINATION))
        {
            positions.add(node.getPosition());
        }
        for(int i = 0; i < node.getChildren().size(); i++)
        {
            positions.addAll(getPositions(node.getChildren().get(i)));
        }
        return positions;
    }

    private String markPatternInText(Integer startPosition, String pattern)
    {
        String matchingTextLHS = mFullText.substring(0, startPosition);
        String matchingText = mFullText.substring(startPosition, startPosition + pattern.length());
        String matchingTextRHS = mFullText.substring(startPosition + pattern.length());
        return matchingTextLHS + "[" + matchingText + "]" + matchingTextRHS;
    }

    private void addChildNode(Node parentNode, String text, int position)
    {
        parentNode.getChildren().add(new Node(text, position, parentNode));
    }

    private void extendNode(Node node, String newText, int position)
    {
        String currentText = node.getText();
        String commonPrefix = getLongestCommonPrefix(currentText, newText);

        if(!commonPrefix.equals(currentText))
        {
            String parentText = currentText.substring(0, commonPrefix.length());
            String childText = currentText.substring(commonPrefix.length());
            splitNodeToParentAndChild(node, parentText, childText);
        }

        String remainingText = newText.substring(commonPrefix.length());
        addChildNode(node, remainingText, position);
    }

    private void splitNodeToParentAndChild(Node parentNode, String parentNewText, String childNewText)
    {
        Node childNode = new Node(childNewText, parentNode.getPosition(), parentNode);

        if(parentNode.getChildren().size() > 0)
        {
            while(parentNode.getChildren().size() > 0)
            {
                childNode.getChildren().add(parentNode.getChildren().remove(0));
            }
        }

        parentNode.getChildren().add(childNode);
        parentNode.setText(parentNewText);
        parentNode.setPosition(POSITION_UNDEFINED);
    }

    private String getLongestCommonPrefix(String str1, String str2)
    {
        int compareLength = Math.min(str1.length(), str2.length());
        for(int i = 0; i < compareLength; i++)
        {
            if(str1.charAt(i) != str2.charAt(i))
            {
                return str1.substring(0, i);
            }
        }
        return str1.substring(0, compareLength);
    }

    private List<Node> getAllNodesInTraversePath(String pattern, Node startNode, boolean isAllowPartialMatch)
    {
        List<Node> nodes = new ArrayList<>();
        for(int i = 0; i < startNode.getChildren().size(); i++)
        {
            Node currentNode = startNode.getChildren().get(i);
            String nodeText = currentNode.getText();
            if(pattern.charAt(0) == nodeText.charAt(0))
            {
                if(isAllowPartialMatch && pattern.length() <= nodeText.length())
                {
                    nodes.add(currentNode);
                    return nodes;
                }

                int compareLength = Math.min(nodeText.length(), pattern.length());
                for(int j = 1; j < compareLength; j++)
                {
                    if(pattern.charAt(j) != nodeText.charAt(j))
                    {
                        if(isAllowPartialMatch)
                        {
                            nodes.add(currentNode);
                        }
                        return nodes;
                    }
                }

                nodes.add(currentNode);
                if(pattern.length() > compareLength)
                {
                    List<Node> nodes2 = getAllNodesInTraversePath(pattern.substring(compareLength), currentNode, isAllowPartialMatch);
                    if(nodes2.size() > 0)
                    {
                        nodes.addAll(nodes2);
                    }
                    else if(!isAllowPartialMatch)
                    {
                        nodes.add(null);
                    }
                }
                return nodes;
            }
        }
        return nodes;
    }

    public String printTree()
    {
        return mRoot.printTree("");
    }

    static class Node
    {
        private String text;
        private int position;
        private final Node parent;
        private final List<Node> children;

        public Node(String word, int position, Node parent)
        {
            this.text = word;
            this.parent = parent;
            this.position = position;
            this.children = new ArrayList<>();
        }

        public String getText()
        {
            return text;
        }

        public void setText(String text)
        {
            this.text = text;
        }

        public int getPosition()
        {
            return position;
        }

        public void setPosition(int position)
        {
            this.position = position;
        }

        public List<Node> getChildren()
        {
            return children;
        }

        public String printTree(String depthIndicator)
        {
            StringBuilder str = new StringBuilder();
            String positionStr = position > -1 ? "[" + position + "]" : "";
            str.append(depthIndicator).append(text).append(positionStr).append("\n");

            for(Node child : children)
            {
                str.append(child.printTree(depthIndicator + "\t"));
            }
            return str.toString();
        }

        @Override
        public String toString()
        {
            return printTree("");
        }
    }
}