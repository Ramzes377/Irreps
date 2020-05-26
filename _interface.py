import partitions
from representations import Representation, SymmetricGroup, AlternatingGroup, get_matrix, cycle_format

from tkinter import Tk
from tkinter import Scrollbar, Menu, Label, Frame,  Text, Entry, Button, Checkbutton, Radiobutton
from tkinter import  END, INSERT, IntVar
from functools import partial
from sympy.combinatorics import Permutation

import sys
import re

STANDART_MESSAGE = 'Please input shape or choose from list after enter n!\nN must be positive integer number!\nShape must be sequence of non increasing integer numbers.\nExamples:\n2, 1\n3, 2\n3, 3, 3\netc...'

try:
    class textContainer(Frame): #Frame 'overload'
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs, borderwidth=1, relief="sunken")

            self.text = Text(self, width=60, height=None, wrap="none", borderwidth=0)

            textVsb = Scrollbar(self, orient="vertical", command=self.text.yview)
            textHsb = Scrollbar(self, orient="horizontal", command=self.text.xview)

            self.text.configure(yscrollcommand=textVsb.set, xscrollcommand=textHsb.set)
            self.text.grid(row=0, column=0, sticky="nsew")

            textVsb.grid(row=0, column=1, sticky="ns")
            textHsb.grid(row=1, column=0, sticky="ew")

            self.text.insert(INSERT, STANDART_MESSAGE)

            self.grid_rowconfigure(0, weight=1)
            self.grid_columnconfigure(0, weight=1)


        def insert(self, *args, **kwargs):
            self.text.insert(*args, **kwargs)

        def delete(self, *args, **kwargs):
            self.text.delete(*args, **kwargs)

    def get_repres_instance():
        group_shape = None
        group = SymmetricGroup if switch_variable.get() else AlternatingGroup
        textContainer.delete('1.0', END)
        try:
            temp = shape_widget.get()
            group_shape = tuple(map(int, temp.split(',')))
            n_widget.delete(0, END)
            n_widget.insert(0, sum(group_shape))
        except:
            textContainer.insert(INSERT, "ERROR! Incorrect shape input!\n\n\n" + STANDART_MESSAGE)

        repres = Representation(_group=group, shape = group_shape, get_exact_output= exact_output.get())
        repres.find_simple_transps_reprs()
        return repres

    def calculate():
        repres = get_repres_instance()

        if elems_find.get():
            repres.find_remain_elems()

        if file_out.get(): #write in file
            out, filename = repres.write_to_file(console_output = console_out.get())
            textContainer.insert(INSERT, out)
            textContainer.insert(INSERT, f'\nSuccesfully writed in file {filename}.txt')
        elif console_out.get(): #write in console
            textContainer.insert(INSERT, repres.output)

    def on_n_input(event):
        def throw_shape(shape):
            shape_widget.delete(0, END)
            shape_widget.insert(0, shape)

        try:
            group_n = int(n_widget.get())
        except ValueError:
            n_widget.delete(0, END)
            n_widget.insert(0, 'Invalid input')
            return
        menu.delete(0, END)
        for p in partitions.partitions(group_n, True):
            label = ','.join(map(str, p.sequence))
            label += '*' if p.is_self_conjugated else ' '
            f = partial(throw_shape, label[:-1])
            menu.add_command(label=label, command= f)
        menu.post(n_widget.winfo_rootx(), n_widget.winfo_rooty())

    def on_perm_input(event):
        repres = get_repres_instance()

        n = int(n_widget.get())

        str_perm = perm_widget.get()
        subperm = re.split(r'[)(]', str_perm)
        t = [tuple(map(int, s.split(','))) for s in subperm if s]

        perm = Permutation(t, size=n + 1)

        # perm_widget.delete(0, END)
        # perm_widget.insert(0, str(perm))

        textContainer.delete('1.0', END)

        matrix = repres[perm]
        out = cycle_format(perm) + get_matrix(matrix, exact_output.get())

        if file_out.get(): #write in file
            filename = repr(repres.group) + '-' + str(repres.group._shape) + '-' + cycle_format(perm)[1:] + '.txt'
            with open(filename, 'w', encoding='utf8') as f:
                f.write(out)
            if console_out.get():
                textContainer.insert(INSERT, out)
            textContainer.insert(INSERT, f'\nSuccesfully writed in file {filename}')
        elif console_out.get(): #write in console
            textContainer.insert(INSERT, out)

    def window_close():
        root.quit()
        root.destroy()
        sys.exit()

except Exception as e:
    textContainer.delete('1.0', END)
    error = f'Oh looks like we got some issue:\n{e}'
    textContainer.insert(INSERT, error)

if __name__  == "__main__":
    root = Tk()
    root.title('Representations')
    root.protocol('WM_DELETE_WINDOW', window_close)
    root.resizable(False, False)

    menu = Menu(root, tearoff=0) #menu-container for partitions choice

    switch_variable = IntVar(value=1) #group choice
    sym_group_button = Radiobutton(root, text="Symmetric Group",
                                   variable=switch_variable,
                                   indicatoron=False, value=1,
                                   width=20,
                                   background = '#d1d1a1').grid(row = 1, column = 0)
    alt_group_button = Radiobutton(root, text="Alternating Group",
                                   variable=switch_variable,
                                   indicatoron=False,
                                   value=0,
                                   width=20,
                                   background = '#cde391').grid(row = 1, column = 1)

    n_label = Label(root, text="N: ").grid(row=2, column=0) #n of group field
    n_widget = Entry(root)
    n_widget.grid(row=2, column=1)
    n_widget.bind("<Return>", on_n_input)

    shape_label = Label(root, text="Shape (partition): ").grid(row=3, column=0) #shape of group
    shape_widget = Entry(root)
    shape_widget.grid(row=3, column=1)

    console_out = IntVar(0) #drop output in textContainer
    write = Checkbutton(root, text="Console output", variable = console_out, onvalue = 1, offvalue = 0)
    write.grid(row = 4, column = 0)

    file_out = IntVar(0) #drop output in file
    in_file = Checkbutton(root, text="Write in file", variable = file_out, onvalue = 1, offvalue = 0)
    in_file.grid(row = 4, column = 1)

    elems_find = IntVar(0) #find only generators or all elements representations
    only_gens = Radiobutton(root, text="Calculate only generators\n(simple transpositions)", variable=elems_find,
                                indicatoron=True, value=0, width=20).grid(row = 5, column = 0)
    all_elems = Radiobutton(root, text="Calculate all elements\n(may take long time)", variable=elems_find,
                                 indicatoron=True, value=1, width=20).grid(row = 5, column = 1)

    exact_output = IntVar(0) #find representations in "exact" mode or not
    float_out = Radiobutton(root, text="float output", variable=exact_output,
                                indicatoron=True, value=0, width=20).grid(row = 6, column = 0)
    exact_out = Radiobutton(root, text="exact output\n(support square roots)", variable=exact_output,
                                 indicatoron=True, value=1, width=20).grid(row = 6, column = 1)

    perm_label = Label(root, text="Calculate permutation: ").grid(row=7, column=0) #permutation input field
    perm_widget = Entry(root)
    perm_widget.grid(row=7, column=1)
    perm_widget.bind("<Return>", on_perm_input)

    b = Button(root, text = 'Calculate!', command = calculate, bg = '#639dc2', width = 40).grid(row = 8, column = 0, columnspan=2)

    textContainer = textContainer(root)
    textContainer.grid(row=0, column=2, rowspan=9, columnspan=10)

    root.mainloop()
