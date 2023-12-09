import markdown


class MarkdownBuilder(str):
    def __init__(self) -> None:
        self.content = ""
    def add_header(self, header, level = '#'):
        self.content += f"{level} {header}\n\n"
    
    def add_codeblock(self):
        self.content += "```"
        self.content += "\n"
        
    def close_codeblock(self):
        self.content += "```"
        self.content += "\n\n"
        
    def add_codeline(self, line):
        self.content += ">>>" + " " + line + "\n"
    
    def add_text(self, text):
        self.content += text + "\n\n"
        
    def add_empty_line(self):
        self.content += "\n"
        
    def add_code_intend(self, code):
        self.content += "..." +"\t\t" +  code + "\n"
        

    
    