print("Hey I'm Glitch, your new assistant")
x= input("What's your name:")
print ("Hey there ", x)
y=input("What's your gender?")
if y== "male":
    print("welcome gentleman")
elif y== "female":
    print("Hey sweet lady/girl")
else:
    print("so you don't prefer to say")
z=input("Where are you from :")
if z=="hyderabad":
    print("I'm sure that you've visited charminar atleast once")
elif z=="manuguru":
    print("I've heard that there is a huge traffic of lizards, isn't it?")
elif z=="mumbai":
    print("so you might have eaten the most delicious vada pav")
elif z=="delhi":
    print("I pity your lungs")
else:
    print("I'm sure that",z,"is a safe place")
a= input("are you a math student or a science student")
if a=="math":
    print("Oh you are a mathie")
    b=input( "2+ 3.5 =?")
    if b== "5.5":
        print("You are too good at math")
    else:
        print("You are poor at math")
    f = input("Do you want me to ask some more questions?:")
    if f == "yes":
        print("Um sure")
        w = (input("what is log2 (upto 4 digits) "))
        if w == "0.3010":
            print("Superb")
        else:
            print("You need to improve")
    else:
         print("OK fine")
elif a=="science":
    print("Oh you gonna be a doctor or an engineer")
    c=input("What is the unit of magnetic flux density?")
    if c== "tesla":
        print("Excellent !!")
    else:
        print("You need to improve at physics")
    d= input("To which chamber of heart does pure blood enter first?")
    if d== "right auricle" or "right atrium":
        print("Good luck future doctor!!")
    else:
        print("concentrate on more on biology")

print ("Its very nice  talking to you")
e= input("Did you enjoy the conversation with glitch?:")
if e== "yes":
    print("That's very kind of you")
else:
    print("Sorry, we'll try to bring up more fun stuff")



