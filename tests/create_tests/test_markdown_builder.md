


## Launch Markdown Builder

```
>>> from src.functional_testing.create_tests.markdown_builder import MarkdownBuilder
>>> markdown_builder = MarkdownBuilder()

```

### add header

```
>>> markdown_builder.add_header('header')
>>> assert markdown_builder.content == "# header\n\n"

```

### add codeblock

```
>>> markdown_builder.content = ""
>>> markdown_builder.add_codeblock()
>>> assert markdown_builder.content == "```\n"

```

### close_codeblock

```
>>> markdown_builder.content = ""
>>> markdown_builder.close_codeblock()
>>> assert markdown_builder.content == "```\n\n"

```

### add codeline

```
>>> markdown_builder.content = ""
>>> markdown_builder.add_codeline('codeline')
>>> assert markdown_builder.content == ">>> codeline\n"

```

### add text

```
>>> markdown_builder.content = ""
>>> markdown_builder.add_text('text')
>>> assert markdown_builder.content == "text\n\n"

```

### add empty line

```
>>> markdown_builder.content = ""
>>> markdown_builder.add_empty_line()
>>> assert markdown_builder.content == "\n"

```
